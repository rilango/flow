import logging
from .base import BaseTask

log = logging.getLogger(__name__)


class Condition(BaseTask):

    def __call__(self, func, **kwargs):

        task_args = {}
        task_args.update(vars(self))
        task_args.update(kwargs)

        return ExecutableCondition(func, **task_args)


class ExecutableCondition(Condition):

    def __init__(self, func, **kwargs):
        Condition.__init__(self, **kwargs)
        assert callable(func), f'Condition "{self.name}" is not callable'
        self.func = func

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)
