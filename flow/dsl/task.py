import logging

from flow.config import Configuration

from .base import BaseTask

log = logging.getLogger(__name__)


class Task(BaseTask):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        assert self.inputs is not None, f'Task "{self.name}" must have inputs'
        assert self.outputs is not None, f'Task "{self.name}" must have outputs'

        # Just to make life easier
        if self.outputs is not None and type(self.outputs) is not list:
            self.outputs = [self.outputs]

    def __call__(self, func, **kwargs):
        BaseTask.__call__(self, **kwargs)
        task_args = self.copy_attributes(**kwargs)
        return ExecutableTask(func, **task_args)


class ExecutableTask(BaseTask):

    def __init__(self, func, **kwargs):
        BaseTask.__init__(self, **kwargs)
        assert callable(func), f'Task "{self.name}" is not callable'
        self.func = func

    async def __call__(self, **kwargs):
        inputs = {k: kwargs[k] for k in self.inputs if kwargs.get(k) is not None}
        assert len(self.inputs) == len(inputs), \
            f'Invalid # of inputs to {self.name} expected {", ".join(inputs)}'

        log.debug(f'Running task "{self.name}({self.id})"')
        log.debug(f'Task arguments {kwargs}')

        # Update config with task details. This includes everything from the
        # Tasks, ExecuteTask, etc..
        task_args = self.copy_attributes(**kwargs)
        del task_args['func']
        task_args['task_id'] = self.id

        # TODO: Revisit to return result with some config. The config must include
        #       - task_id
        try:
            return self.func(config=Configuration(), **task_args)
        except Exception as e:
            log.error(f'Task "{self.name}({self.id})" failed: {e}')
            raise e
