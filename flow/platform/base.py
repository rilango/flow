import logging

log = logging.getLogger(__name__)


class BasePlatform(object):
    '''
    Base class for all platform implementations
    '''

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for starting a defining a platform
        '''
        self.ws_root = kwargs.get('platform_ws_root')

        assert self.ws_root is not None, \
            f'Platform needs ws_root.'

    def create_ws(self, **kwargs):
        '''
        Creates a workspace for an instance of pipeline on the platform.

        Returns:
            str: Path/resource location for the job to store artifacts.
        '''
        return NotImplemented

    def create_sub_ws(self, parent_ws, **kwargs):
        '''
        Creates a workspace for a task with-in a job.

        Returns:
            str: Path/resource location for the task to store artifacts.
        '''
        return NotImplemented


    def morph_arguments(self, config, **kwargs):
        '''
        Morphs the arguments for a task.

        Args:
            config (dict): config of the task.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            dict: Config of the task.
        '''
        return NotImplemented

class BasePlatformTask(object):
    '''
    Base class for all Task implementations
    '''

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for starting a defining a platform
        '''
        self.params = kwargs.get('params')
        self.task_workspace = kwargs.get('task_workspace')

        assert self.inputs is not None, f'Task "{self.name}" must have inputs'
        assert self.outputs is not None, f'Task "{self.name}" must have outputs'

    def run(self):
        '''
        Runs the containerized task
        '''
        pass
