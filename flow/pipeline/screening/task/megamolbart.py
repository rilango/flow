import logging

from flow.platform import ContainerizedTask
from flow.platform.container import LoadBalancedService
from flow.utils.megamolbart import sample

log = logging.getLogger(__name__)


class MegaMolBARTService(ContainerizedTask):

    container_image = 'nvcr.io/nv-drug-discovery-dev/megamolbart:latest'
    entry_point =\
        '''\
        bash -c "cd /opt/nvidia/cheminfomatics/megamolbart && python3 -m megamolbart"
        '''

    def __init__(self, **kwargs):
        self._service_port = kwargs.get('service_port', '127.0.0.1:50051')
        ContainerizedTask.__init__(self, **kwargs)

    def __call__(self, **kwargs):
        log.info(f'Running Molbart service: {self.container_image}')
        ContainerizedTask.__call__(self, **kwargs)

    def is_available(self):
        try:
            sample('CC(=O)O',
                   radius=1,
                   num_sample=1,
                   service_port=self._service_port)
            return True
        except Exception as e:
            log.exception(e)
            return False


class LBMegaMolBARTService(LoadBalancedService):

    work_folder = './support/docker/megamolbart'
    compose_config = 'support/docker/megamolbart/docker-compose.yml'
    lb_config = './support/docker/megamolbart/nginx/nginx.conf'

    def __init__(self, **kwargs):
        #TODO: Service port is hardcoded for now due to a bug in docker-compose
        #      https://github.com/docker/compose/issues/7928
        self._service_port = kwargs.get('service_port', 'localhost:50052')
        LoadBalancedService.__init__(self, **kwargs)

    def __call__(self, **kwargs):
        log.info(f'Running Molbart load balanced service: {self.compose_config}')

        LoadBalancedService.__call__(self,
                                     service_port=self._service_port,
                                     nginx_config=self.lb_config,
                                     scale_service={'megamolbart': 4},
                                     **kwargs)

    def is_available(self):
        try:
            log.info(f'Checking if Molbart service is available {self._service_port}')
            sample('CC(=O)O',
                   radius=1,
                   num_sample=1,
                   service_port=self._service_port)
            return True
        except Exception as e:
            log.exception(e)
            return False
