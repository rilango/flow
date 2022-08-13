import logging

from flow.pipeline.screening.data import ScreeningData
from flow.platform import ContainerizedTask


log = logging.getLogger(__name__)


class SimpleScreen():

    input_spec = ['generation']
    output_spec = ['curr_generation']

    def __call__(self, **kwargs):

        generation = kwargs['generation']
        num_mol = kwargs['top_k']

        return ScreeningData().promote_smiles(generation, num_mol)


class OptimizedScreen(ContainerizedTask):

    input_spec = ['generation', 'pdbqt_file', 'score_config']
    output_spec = ['z']
    container_image = 'mgltools:latest'
    entry_point = '''\
        /opt/conda/envs/rapids/bin/python\
            /code/flow/pipeline/screening/optimization.py\
                -g {{generation}}\
                -n {{top_k}}\
                --db /pipeline/common.sqlite3
        '''

    def __call__(self, **kwargs):
        generation = kwargs['generation']
        log.info(f'Screening molecules from generation {generation}')

        ContainerizedTask.__call__(self, **kwargs)

        return generation
