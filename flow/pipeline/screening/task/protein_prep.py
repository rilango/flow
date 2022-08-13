import logging
from flow.platform import ContainerizedTask
from flow.pipeline.screening.data import ScreeningData


log = logging.getLogger(__name__)


class PrepareProtein(ContainerizedTask):

    input_spec = ['protein_file', 'pp_repair_type']
    output_spec = ['pdbqt_file']
    container_image = 'mgltools:latest'
    entry_point =\
        '''\
            bash -c \
            "
            conda run -n mgltools \
                python2 /opt/conda/envs/mgltools/bin/prepare_receptor4.py \
                    -r {{protein_file}} \
                    -A {{pp_repair_type}} \
                    -o {{protein_file}}qt"
        '''

    def __call__(self, **kwargs):
        data = ScreeningData()
        data.initialize_db()
        protein_file = kwargs['protein_file']

        log.info(f'Creating pdbqt file for {protein_file}')
        ContainerizedTask.__call__(self, **kwargs)
        return f'{protein_file}qt'