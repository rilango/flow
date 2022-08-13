import io
import os
import pytest
import logging
import asyncio

from flow.dsl import task, service, pipeline, loop, condition, when
from flow.config import Configuration
from flow.pipeline.screening.task import (PrepareProtein,
                                         MegaMolBARTService,
                                         GenerateMolecule,
                                         ScoreMolecule,
                                         ScreenMolecule)


log = logging.getLogger(__name__)


def iterate_five(**kwargs):
    return kwargs['_loop_count'] < 20


def next_generation(**kwargs):
    kwargs['generation'] = kwargs['generation'] + 1
    return kwargs


@pytest.mark.skip(reason="This is used for verification and benchmarking only")
def test_vs_pipeline():

    config = '''
    platform_type: containers
    platform_ws_root: /raid/drugdiscovery/cheminformatics/vs_screening
    platform_db: /raid/drugdiscovery/cheminformatics/vs_screening/{{pipeline_id}}/common.sqlite3
    platform_db_link: /pipeline/common.sqlite3
    '''
    data = io.StringIO()
    data.write(config)
    config = Configuration()
    config.load_config(data)
    log.info(f'Starting pipeline with config: {config.config}')

    mol_generator = GenerateMolecule()
    virtual_screening =\
        pipeline(name='virtual_screening')(
            service(name='mega_molbart')(MegaMolBARTService(ports=['50051:50051'])),
            task(name='prepare_protein',
                 inputs=['protein_file', 'pp_repair_type', 'score_config'],
                 outputs=['pdbqt_file'])(PrepareProtein()),
            task(name='GenerateMolecule',
                 inputs=['smiles', 'generation', 'min_lipinsky_score'],
                 outputs=['curr_generation'],
                 constants={'service_port': '127.0.0.1:50051',
                         'num_sample': 20})(mol_generator),
            loop(name='Check Pose Score', increment=next_generation)(
                task(name='ScoreMolecule',
                     inputs=['generation', 'pdbqt_file'],
                     outputs=['status'])(ScoreMolecule()),
                when(name='GenerateWhen')(
                    condition(name='GenerateWhen',
                              inputs='_loop_count')(iterate_five),
                    task(name='ScreenMolecule',
                         inputs=['generation', 'top_k'],
                         outputs=['generation'])(ScreenMolecule()),
                    task(name='GenerateMolecule',
                         inputs=['generation', 'min_lipinsky_score'],
                         outputs=['test'],
                         constants={'service_port': '127.0.0.1:50051',
                                    'num_sample': 20,
                                    'use_generation': True})(mol_generator),
                ),
                condition(name='Check Pose Score',
                          inputs='_loop_count')(iterate_five),
            )
        )

    pipeline_input = {
                      'protein_file': 'file://./tests/data/6y2g/rec.pdb',
                      'score_config': 'file://./tests/data/6y2g/config',
                      'pp_repair_type': 'checkhydrogens',
                      'smiles': 'CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@@](C#N)(c2ccc3c(N)ncnn23)[C@H](O)[C@@H]1O)Oc1ccccc1',
                      'generation': 0,
                      'top_k': 3,
                      'min_lipinsky_score': 0}
    asyncio.run(virtual_screening(**pipeline_input))

    ws_root = config.get('platform_ws_root')

    ws = os.path.join(ws_root, virtual_screening.id)
    assert os.path.exists(ws), 'Workspace was not created'
    log.info(f'Pipeline {virtual_screening.id} workspace was at {ws}')

    assert os.path.exists(os.path.join(ws, '_spec.json')), \
        'Workspace spec file was not created'

    with open(os.path.join(ws, '_spec.json'), 'r') as f:
        log.debug(f.read())
