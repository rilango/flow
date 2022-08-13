#!/usr/bin/env python3
'''
For a given list of smiles, generate conformers and store them in the database.

This scripts works in teandum with flow.pipeline.screening package.
'''
import re
import os
import sys
import uuid
import argparse
import logging
import asyncio
import concurrent
import multiprocessing

from functools import partial
from subprocess import run

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdmolfiles

from flow.pipeline.screening.data import ScreeningData


logging.basicConfig(format='%(asctime)s %(name)s [%(levelname)s]: %(message)s',
                    level=logging.WARNING)
log = logging.getLogger(__name__)

def argparser(argv):
    parser = argparse.ArgumentParser(description='Parse ligand preparation arguments')

    parser.add_argument('-g', '--generation',
                        dest='generation',
                        type=int,
                        required=True,
                        help='Generation of the ligand library')

    parser.add_argument('-r', '--receptor',
                        dest='receptor',
                        type=str,
                        required=True,
                        help='Receptor file')

    parser.add_argument('-c', '--score_config',
                        dest='score_config',
                        type=str,
                        required=True,
                        help='Score config with binding site and box size for vina')

    parser.add_argument('--cpus',
                        dest='cpus',
                        type=int,
                        default=None,
                        required=False,
                        help='Number of CPUs to use')

    parser.add_argument('--db',
                        dest='db_url',
                        type=str,
                        default=None,
                        required=False,
                        help='Database connection string.')

    args = parser.parse_args(argv)
    return args


def score_pose(score_file):
    with open(score_file, 'r') as fh:
        lines = fh.read()
        scorelines = re.findall(r'REMARK VINA RESULT.*', lines)
        min_score = sys.maxsize
        score_model = None
        cnt = 1
        for scoreline in scorelines:
            score = float(scoreline.split()[3])
            if min_score > score:
                min_score = score
                score_model = cnt
            cnt += 1

        return min_score, score_model


def generate_pose(receptor_file, score_config, ligand_file, gsmiles_id, cid):
    ligand_folder = os.path.dirname(ligand_file)

    cpu_cnt = os.cpu_count()//2
    out_file = f'{ligand_folder}/{gsmiles_id}_{cid}_vina.pdbqt'
    log_file = f'{ligand_folder}/{gsmiles_id}_{cid}.log'
    cmd = ['vina',
           '--receptor', receptor_file, \
           '--ligand', ligand_file, \
           '--out', out_file, \
           '--log', log_file, \
           '--cpu', str(cpu_cnt), \
           '--config', score_config]
    log.info(' '.join(cmd))
    result = run(' '.join(cmd), capture_output=True, shell=True)
    if result.returncode == 0:
        return score_pose(out_file)
    else:
        log.error(f'Failed to generate pose for {gsmiles_id}_{cid} with error {result.stderr}')
        return (None, None)


def score_conformer(receptor_file, score_config, gsmiles_id, mol, conformer_file, cid):
    rdmolfiles.MolToPDBFile(mol, conformer_file, confId=cid)
    data_dir = os.path.dirname(conformer_file)
    cmd = f'''
        cd {data_dir};
        python2 /opt/conda/envs/mgltools/bin/prepare_ligand4.py \
            -l {conformer_file} \
            -o {conformer_file}qt \
            -A bonds_hydrogens
        '''
    result = run(cmd, capture_output=True, shell=True)
    if result.returncode == 0:
        min_score, score_model = generate_pose(receptor_file,
                                               score_config,
                                               f'{conformer_file}qt',
                                               gsmiles_id,
                                               cid)
        return (min_score, score_model)
    else:
        log.error(f'Failed to prepare ligand for {gsmiles_id}_{cid} with error {result.stderr}')
        return (None, None)


def _score_conformer(receptor_file,
                     score_config,
                     gsmiles_id,
                     mol,
                     conformer_file,
                     cid,
                     db_file):
    '''
    Score conformers and store them in the database.
    '''
    min_score, score_model = score_conformer(receptor_file,
                                             score_config,
                                             gsmiles_id,
                                             mol,
                                             conformer_file,
                                             cid)
    if min_score is not None:
        log.info(f'Score for {gsmiles_id}_{cid} is {min_score} from model {score_model}')
        ScreeningData(datasource=db_file)\
            .add_conformer_score(gsmiles_id, cid, min_score, score_model)


def generate_conformers(gsmiles_id, smiles):
    log.info(f'Computing conformers for {smiles}')
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        log.info(f'Computing conformers failed for {smiles}')
        return
    mol = Chem.AddHs(mol, addCoords=True)

    num_conformers=10

    params = rdDistGeom.ETKDGv2()
    params.pruneRmsThresh = 0.01
    params.randomSeed = 42
    params.numThreads = 0
    params.clearConfs = True
    params.maxIterations = 1000

    conformers = rdDistGeom.EmbedMultipleConfs(mol,
                                               numConfs=num_conformers,
                                               params=params)
    return {'conformers': conformers,
            'gsmiles_id': gsmiles_id,
            'mol': mol}


async def _produce_conformers(queue, generation, db_file):
    all_smiles = ScreeningData(datasource=db_file).\
        fetch_gsmiles_by_generation(generation,
                                    include_unscored=True)
    log.info(f'Number of molecules in gen {generation} is {len(all_smiles)}')

    for gsmiles_id, details in all_smiles.items():
        conformer_spec = generate_conformers(gsmiles_id, details['smiles'])
        log.info(f'Adding conformer to queue')
        await queue.put(conformer_spec)


async def process_conformers(queue, receptor_file, score_config, db_file, cpu_cnt):
    conformer_dir = os.path.join('/pipeline', 'ligand')
    os.makedirs(conformer_dir, exist_ok=True)

    if cpu_cnt is None or cpu_cnt == 0:
        cpu_cnt = multiprocessing.cpu_count() // 2

    loop = asyncio.get_running_loop()
    with concurrent.futures.ThreadPoolExecutor(cpu_cnt) as pool:
        while True:
            log.info(f'Waiting for conformer spec...')
            conformer_spec = await queue.get()
            gsmiles_id = conformer_spec['gsmiles_id']
            mol = conformer_spec['mol']
            log.info(f'Processing conformer for gsmiles_id {gsmiles_id}...')
            for cid in conformer_spec['conformers']:
                conformer_file = f'{conformer_dir}/{gsmiles_id}_{cid}.pdb'

                partial_func = partial(_score_conformer,
                                       receptor_file,
                                       score_config,
                                       gsmiles_id,
                                       mol,
                                       conformer_file,
                                       cid,
                                       db_file)
                results = loop.run_in_executor(pool, partial_func)
            queue.task_done()

def score_molecule(smiles,
                   receptor_file,
                   score_config,
                   conformer_dir='/tmp'):
    os.makedirs(conformer_dir, exist_ok=True)

    gsmiles_id = uuid.uuid1().int
    conformer_spec = generate_conformers(gsmiles_id, smiles)
    log.debug(f'Adding conformer to queue{conformer_spec}')
    gsmiles_id = conformer_spec['gsmiles_id']
    mol = conformer_spec['mol']

    log.info(f'Processing conformer for gsmiles_id {gsmiles_id}...')

    cpu_cnt = 4
    if cpu_cnt is None or cpu_cnt == 0:
        cpu_cnt = multiprocessing.cpu_count() // 2
    partial_func = partial(score_conformer,
                           receptor_file,
                           score_config,
                           gsmiles_id,
                           mol)

    with concurrent.futures.ThreadPoolExecutor(max_workers=cpu_cnt) as executor:
        futures = {executor.submit(partial_func, f'{conformer_dir}/{gsmiles_id}_{cid}.pdb', cid): \
            cid for cid in conformer_spec['conformers']}

        min_score = None
        score_model = None
        for future in concurrent.futures.as_completed(futures):
            try:
                score, score_model = future.result()
                if min_score is None or score < min_score:
                    min_score = score
                    score_model = score_model
            except Exception as exc:
                log.warning(f'{future} generated an exception: {exc}')

    return min_score, score_model

async def main():
    args = argparser(sys.argv[1:])
    db_url = args.db_url
    generation = args.generation
    receptor_file = args.receptor
    score_config = args.score_config
    cpus = args.cpus

    queue = asyncio.Queue()
    consumer = asyncio.create_task(
        process_conformers(queue, receptor_file, score_config, db_url, cpus))

    await _produce_conformers(queue, generation, db_url)
    await queue.join()
    consumer.cancel()

if __name__ == '__main__':
    asyncio.run(main())
