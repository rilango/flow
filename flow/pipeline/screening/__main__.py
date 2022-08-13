import sys
import time
import logging
import argparse
import asyncio

from flow.dsl import task, service, pipeline, loop, condition, when
from flow.config import Configuration
from flow.pipeline.screening.data import ScreeningData
from flow.pipeline.screening.task import (PrepareProtein,
                                         MegaMolBARTService,
                                         GenerateMolecule,
                                         ScoreMolecule,
                                         SimpleScreen,
                                         OptimizedScreen)
from flow.pipeline.screening.task.megamolbart import LBMegaMolBARTService


# Log configuration
logging.basicConfig(format='%(asctime)s %(name)s [%(levelname)s]: %(message)s',
                    level=logging.DEBUG,
                    stream = sys.stdout)
log = logging.getLogger(__name__)


def process_args():
    '''
    Process command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Analyze')
    parser.add_argument('--config', '-c',
                        type=str,
                        dest='config',
                        required=True,
                        help='Pipeline configuration file')

    args = parser.parse_args(sys.argv[1:])

    config = Configuration()
    config.load_config(args.config)
    return config


def add_initial_smiles(**kwargs):
    data = ScreeningData()
    data.initialize_db()
    all_smiles = data.add_smiles(kwargs['smiles'], 0)
    return 0


def todo_or_nottodo(**kwargs):
    return kwargs['_loop_count'] < 10


def next_generation(**kwargs):
    kwargs['generation'] = kwargs['generation'] + 1
    return kwargs


def simple_pipeline():
    mol_generator = GenerateMolecule()
    virtual_screening =\
        pipeline(name='virtual_screening')(
            service(name='mega_molbart')\
                (MegaMolBARTService(ports=['50051:50051'], service_port='127.0.0.1:50051')),
            task(name='prepare_protein',
                 inputs=['protein_file', 'pp_repair_type', 'score_config'],
                 outputs=['pdbqt_file'])(PrepareProtein()),
            task(name='Add initial smiles',
                 inputs=['smiles'],
                 outputs=['generation'])(add_initial_smiles),
            loop(name='Check Pose Score', increment=next_generation)(
                condition(name='Check Pose Score',
                          inputs='_loop_count')(todo_or_nottodo),
                task(name='GenerateMolecule',
                     inputs=['generation', 'min_lipinsky_score'],
                     outputs=['curr_generation'],
                     constants={'service_port': '127.0.0.1:50051',
                                'num_sample': 100})(mol_generator),
                task(name='ScoreMolecule',
                     inputs=['curr_generation', 'pdbqt_file'],
                     outputs=['status'])(ScoreMolecule()),
                when(name='ScreenWhen')(
                    condition(name='ScreenWhen',
                              inputs='_loop_count')(todo_or_nottodo),
                    task(name='ScreenMolecule',
                         inputs=['curr_generation', 'top_k'],
                         outputs=['generation'])(SimpleScreen()),
                ),
            )
        )

    return virtual_screening


def simple_rbf_pipeline():
    mol_generator = GenerateMolecule()
    virtual_screening =\
        pipeline(name='virtual_screening')(
            service(name='mega_molbart')(
                MegaMolBARTService(ports=['50051:50051'],
                                   service_port='127.0.0.1:50051')),
            task(name='prepare_protein',
                 inputs=['protein_file', 'pp_repair_type', 'score_config'],
                 outputs=['pdbqt_file'])(PrepareProtein()),
            task(name='Add initial smiles',
                 inputs=['smiles'],
                 outputs=['generation'])(add_initial_smiles),
            task(name='GenerateMolecule',
                 inputs=['generation', 'min_lipinsky_score'],
                 outputs=['curr_generation'],
                 constants={'service_port': '127.0.0.1:50051',
                            'num_sample': 10,
                            'span': 10})(mol_generator),
            loop(name='Check Pose Score', increment=next_generation)(
                condition(name='Check Pose Score',
                          inputs='_loop_count')(todo_or_nottodo),
                task(name='ScoreMolecule',
                     inputs=['curr_generation', 'pdbqt_file'],
                     outputs=['status'])(ScoreMolecule()),
                task(name='ScreenMolecule',
                     inputs=['curr_generation', 'top_k'],
                     outputs=['generation'])(OptimizedScreen()),
            )
        )

    return virtual_screening


def simple_regression_pipeline():
    mol_generator = GenerateMolecule()
    virtual_screening =\
        pipeline(name='virtual_screening')(
            service(name='mega_molbart')\
                (MegaMolBARTService(ports=['50051:50051'], service_port='127.0.0.1:50051')),
            task(name='prepare_protein',
                 inputs=['protein_file', 'pp_repair_type', 'score_config'],
                 outputs=['pdbqt_file'])(PrepareProtein()),
            task(name='Add initial smiles',
                 inputs=['smiles'],
                 outputs=['generation'])(add_initial_smiles),
            loop(name='Check Pose Score', increment=next_generation)(
                condition(name='Check Pose Score',
                          inputs='_loop_count')(todo_or_nottodo),
                task(name='GenerateMolecule',
                     inputs=['generation', 'min_lipinsky_score'],
                     outputs=['curr_generation'],
                     constants={'service_port': '127.0.0.1:50051',
                                'num_sample': 20})(mol_generator),
                task(name='ScoreMolecule',
                     inputs=['curr_generation', 'pdbqt_file'],
                     outputs=['status'])(ScoreMolecule()),
                when(name='ScreenWhen')(
                    condition(name='ScreenWhen',
                              inputs='_loop_count')(todo_or_nottodo),
                    task(name='ScreenMolecule',
                         inputs=['curr_generation', 'top_k'],
                         outputs=['generation'])(OptimizedScreen()),
                ),
            )
        )

    return virtual_screening


def main():
    '''
    Entry point to execute pipeline
    '''
    config = process_args()
    log.info(f'Starting pipeline with config: {config.config}')

    # virtual_screening = simple_pipeline()
    virtual_screening = simple_rbf_pipeline()
    # virtual_screening = simple_regression_pipeline()
    # virtual_screening = simple_test()

    pipeline_input = config.get('inputs')

    start = time.time()
    asyncio.run(virtual_screening(**pipeline_input))
    log.info(f'Runtime {time.time() - start}')


if __name__ == '__main__':
    main()
