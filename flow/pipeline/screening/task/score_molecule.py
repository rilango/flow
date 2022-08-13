import logging
from flow.pipeline.screening.data import ScreeningData
from flow.platform import ContainerizedTask


log = logging.getLogger(__name__)


class ScoreMolecule(ContainerizedTask):

    input_spec = ['generation', 'pdbqt_file', 'score_config']
    output_spec = ['z']
    container_image = 'mgltools:latest'
    entry_point = '''\
        /opt/conda/envs/rapids/bin/python\
            /code/flow/pipeline/screening/pose_generate.py\
                -g {{generation}}\
                -r {{pdbqt_file}}\
                -c {{score_config}}\
                --db /pipeline/common.sqlite3
        '''

    def __call__(self, **kwargs):
        generation = kwargs['generation']
        log.info(f'Scoring molecules from generation {generation}')

        ContainerizedTask.__call__(self, **kwargs)

        return ScreeningData().fetch_min_score(generation)
