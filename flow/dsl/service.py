import logging

from .base import BaseTask
from flow.config import Configuration


log = logging.getLogger(__name__)


class Service(BaseTask):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        assert self.inputs is None, f'Task "{self.name}" must not have inputs'
        assert self.outputs is None, f'Task "{self.name}" must not have outputs'
        self._service = True

    def __call__(self, func, **kwargs):
        BaseTask.__call__(self, **kwargs)
        task_args = self.copy_attributes(**kwargs)
        return ExecutableService(func, **task_args)


class ExecutableService(Service):

    def __init__(self, func, **kwargs):
        Service.__init__(self, **kwargs)
        assert callable(func), f'Service "{self.name}" is not callable'
        self.func = func

    def __call__(self, **kwargs):
        assert self.inputs is None and  self.outputs is None, \
            f'Service "{self.name}" cannot have inputs or outputs'
        log.debug(f'Service arguments {kwargs}')

        # Update config with task details. This includes everything from the
        # Tasks, ExecuteTask, etc..
        service_args = self.copy_attributes(**kwargs)
        del service_args['func']
        service_args['service_id'] = self.id

        try:
            self.func(config=Configuration(), **service_args)
        except Exception as e:
            log.error(f'Service "{self.name}({self.id})" failed: {e}')
            raise e
        # return sorted_lists_futures

    def is_available(self):
        if hasattr(self.func, 'is_available'):
            return self.func.is_available()
        else:
            return True
