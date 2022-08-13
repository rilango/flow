import io
import os
import json
import shutil
import logging
import tempfile
from pathlib import Path

from subprocess import PIPE, run
from .base import BasePlatform
from flow.utils.stringutils import substitute

log = logging.getLogger(__name__)


class ContainerPlatform(BasePlatform):

    pipeline_path = '/pipeline'
    task_path = '/task'

    def create_ws(self, **kwargs):
        '''
        Sets up the workspace for a pipeline.

        :param uid: Unique identifier for the pipeline.

        :param platform_ws_root: Root location for the platform
        '''
        uid = kwargs.get('uid')
        ws_root = kwargs.get('platform_ws_root')

        assert uid is not None, 'uid is required.'
        assert ws_root is not None, 'ws_root is required.'

        ws = os.path.join(ws_root, uid)
        os.makedirs(ws)

        with open(os.path.join(ws,  "_spec.json"), "w") as file:
            file.write(json.dumps(kwargs, indent=4))

        return ws

    def create_sub_ws(self, parent_ws, **kwargs):
        pass

    def morph_arguments(self, config, **kwargs):
        '''
        Convert arguments to the target platform.
        - Convert paths to the platform specific path.
        '''
        for k, v in kwargs.items():
            if isinstance(v, io.IOBase) or \
                    (isinstance(v, str) and v.startswith('file://')):
                # input can be a file or a path.
                full_path = v if type(v) == str else v.name
                if full_path.startswith('file://'):
                    full_path = full_path[7:]
                # Copy the file to the platform specific location and replace
                # the arg value with the new path.
                filename = os.path.basename(full_path)
                host_input_dir = os.path.join(
                    config.get_pipeline_path(), 'inputs')

                os.makedirs(host_input_dir, exist_ok=True)
                target_file = os.path.join(config.get_pipeline_path(),
                                           'inputs',
                                           filename)
                shutil.copyfile(full_path, target_file)
                kwargs[k] = os.path.join(
                    self.pipeline_path, 'inputs', filename)

        return kwargs


class ContainerizedTask():
    '''
    Base class for containerized tasks
    '''

    # All containerized tasks must have the following attributes
    input_spec = None
    output_spec = None
    container_image = None
    entry_point = None

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for starting a containerized task.
        '''
        self.volumns = kwargs.get('volumns')
        self.ports = kwargs.get('ports')
        self.params = kwargs.get('params')
        self.work_dir = kwargs.get('work_dir')

    def __call__(self, **kwargs):
        assert kwargs is not None, 'Inputs are required'

        # Set all inputs values
        env_variables = {}
        if kwargs.get('inputs') is not None:
            env_variables = {key: kwargs[key] for key in kwargs['inputs']}

        config = kwargs['config']
        ws_root = config.get('platform_ws_root')
        task_dir = os.path.join(ws_root, *kwargs.get('_node_path'))

        volumns = {config.get_pipeline_path(): config.platform.pipeline_path,
                   task_dir: config.platform.task_path}

        entry_point = substitute(self.entry_point, kwargs)

        # Run the container
        return self.run_container(self.container_image,
                                  entry_point,
                                  env=env_variables,
                                  volumns=volumns)

    def run_container(self, container_image, entrypoint, **kwargs):
        '''
        Execute a container.

        :param container_image: Container image to run.

        :param entrypoint: Entrypoint to run.

        :param env: Environment variables to pass to the container.

        :param task_dir: Working directory to run the container in.

        :param volumns: Volumns to mount to the container.
        '''
        log.debug(
            f'Running container: {container_image} {entrypoint} {kwargs}')

        port_str = ''
        if self.ports is not None:
            for value in self.ports:
                port_str += f' -p {value} '

        env_str = ''
        env = kwargs.get("env")
        if env is not None:
            for key, value in env.items():
                env_str += f' -e {key.upper()}="{value}"'

        mount_str = ''
        volumns = kwargs.get("volumns")
        if volumns is not None:
            for key, value in volumns.items():
                mount_str += f' -v {key}:{value}'

        # Always mount the source directory into the container
        mount_str += f' -v {__file__[: __file__.rfind("flow")]}:/code'
        env_str += f' -e PYTHONPATH="/code:$PYTHONPATH"'

        # TODO: Revisit for resource management
        command = f'''\
            docker run \
                --network=host \
                --gpus all \
                -w /code \
                {env_str} {mount_str} {port_str}\
                {container_image} {entrypoint}
        '''

        log.info(f'Command {command}')
        result = run(command, stdout=PIPE, stderr=PIPE,
                     universal_newlines=True, shell=True)
        log.info(f'Ouput: {result.stdout}\nError: {result.stderr}')
        return result


class LoadBalancedService():
    '''
    Base class for load balanced service
    '''

    work_folder = None
    compose_config = None
    lb_config = None

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for starting a containerized task.
        '''
        pass

    def __call__(self, **kwargs):
        assert kwargs is not None, 'Inputs are required'

        # Run the container
        return self.run_service(**kwargs)

    def run_service(self, **kwargs):
        '''
        Start the service using docker-composer.
        '''
        log.debug(
            f'Starting LB svc: {self.compose_config} {self.lb_config}')

        log.info(f'Creating env file for docker-compose...')
        env_file = tempfile.NamedTemporaryFile(
            prefix=".env.", delete=False)
        for key, value in kwargs.items():
            env_file.write(
                bytes(f'{key.upper()}={str(value)}\n', 'UTF-8'))

        # Service scale
        scale_str = ''
        if kwargs.get('scale_service') is not None:
            for key, value in kwargs['scale_service'].items():
                scale_str += f' --scale {key}={value}'

        # source_root = __file__[: __file__.rfind("flow")]
        # TODO: Revisit for resource management
        command = f'''\
            cd {Path(self.compose_config).parent} &&
            docker-compose --env-file {env_file.name}\
                -f {Path(self.compose_config).name}\
                up {scale_str}
            '''

        log.info(f'Command {command}')
        result = run(command, stdout=PIPE, stderr=PIPE,
                     universal_newlines=True, shell=True)
        log.info(f'Ouput: {result.stdout}\nError: {result.stderr}')
        return result

    def shutdown(self):
        command = f'cd {Path(self.compose_config).parent} && docker-compose down'
        log.info(f'Command {command}')
        result = run(command, stdout=PIPE, stderr=PIPE,
                     universal_newlines=True, shell=True)
        log.info(f'Ouput: {result.stdout}\nError: {result.stderr}')
