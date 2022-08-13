import logging

from .base import TaskBlock
from .pipeline import ExecutablePipeline
from .condition import ExecutableCondition

log = logging.getLogger(__name__)


class When(TaskBlock):
    '''
    Provides a contract to define a Conditional block in a pipeline. This is a
    block that can encapsulate a set of tasks and a condition. The name and
    the condition name to be use for exit critieria must match.
    '''

    def __init__(self, **kwargs):
        TaskBlock.__init__(self, **kwargs)
        self.condition = None

    def __call__(self, *args, **kwargs):
        log.debug(f'Processing "When {self.name}"')
        TaskBlock.__call__(self, *args, **kwargs)

        # Additional validations for loops
        assert self.name in self._graph and \
            isinstance(self._graph[self.name], ExecutableCondition),\
            f'Loop "{self.name}" must contain a matching condition'

        self.inputs, self.outputs = self.effective_values()

        self.condition = self._graph.pop(self.name)
        assert len(self._graph) > 0, 'When must contain at least one task'

        task_args = self.copy_attributes(**kwargs)

        return ExecutableWhen(**task_args)


class ExecutableWhen(ExecutablePipeline):

    def __init__(self, **kwargs):
        ExecutablePipeline.__init__(self, **kwargs)
        self.condition = kwargs.get('condition')
        self.outputs = kwargs.get('outputs')
        self._graph = kwargs.get('_graph')

    async def __call__(self, **kwargs):
        log.debug(f'Running When "{self.name}({self.id})"')
        log.debug(f'When arguments: {kwargs}')
        self._visited_tasks = set()
        self._completed_tasks = set()

        self._state = kwargs.copy()

        # Check if the loop should be terminated
        inputs = self._construct_input(self.condition)

        condition_result = self.condition(**inputs)
        if condition_result:
            await self.execute_block()

        # Filter the output values and return
        return {k: self._state[k] for k in self.outputs}
