import logging

from .base import TaskBlock
from .pipeline import ExecutablePipeline
from .condition import ExecutableCondition

log = logging.getLogger(__name__)


class Loop(TaskBlock):
    '''
    Provides a contract to define a loop in a pipeline. This is a block that can
    encapsulate a set of tasks and a condition. The loop name and the condition
    name to be use for exit critieria must match.
    '''

    def __init__(self, **kwargs):
        TaskBlock.__init__(self, **kwargs)
        self.increment = kwargs.get('increment')
        self.condition = None

    def __call__(self, *args, **kwargs):
        log.debug(f'Processing loop "{self.name}"')
        TaskBlock.__call__(self, *args, **kwargs)

        # Additional validations for loops
        assert self.name in self._graph and \
            isinstance(self._graph[self.name], ExecutableCondition),\
            f'Loop "{self.name}" must contain a matching condition'

        # Add a placeholder node for the loop_count to be replaces at runtime.
        # This will help separate each instance of task in a loop.
        for task in self._graph.values():
            task._node_path.insert(1, f'_loop_{self.id}__')

        self.inputs, self.outputs = self.effective_values()

        self.condition = self._graph.pop(self.name)
        assert len(self._graph) > 0, 'Loop must contain at least one task'

        task_args = self.copy_attributes(**kwargs)
        # increament is a function and needs to be copied as well.
        task_args['increment'] = self.increment

        return ExecutableLoop(**task_args)


class ExecutableLoop(ExecutablePipeline):

    def __init__(self, **kwargs):
        ExecutablePipeline.__init__(self, **kwargs)
        self.condition = kwargs.get('condition')
        self.increment = kwargs.get('increment')
        self.outputs = kwargs.get('outputs')
        self._graph = kwargs.get('_graph')
        self._loop_count = 0

    def _prep_next(self, **kwargs):
        self._visited_tasks = set()
        self._completed_tasks = set()

        # For next iteration of loop merge the state and kwargs
        merged = {}
        merged.update(self._state)
        merged.update(kwargs)
        return merged

    async def __call__(self, **kwargs):
        self._loop_count += 1

        log.debug(f'Running Loop "{self.name}({self.id})" - iter {self._loop_count}')
        log.debug(f'Loop arguments: {kwargs}')

        self._state = kwargs.copy()
        self._state['_loop_count'] = self._loop_count

        # Update the loop path placeholder with the loop count.
        for task in self._graph.values():
            task._node_path[-2] = f'_loop_{self.id}_{self._loop_count}'

        await self.execute_block()

        # Check if the loop should be terminated
        inputs = self._construct_input(self.condition)

        if self.increment:
            kwargs = self.increment(**kwargs)
        condition_result = self.condition(**inputs)
        if condition_result:
            merged_data = self._prep_next(**kwargs)
            await self(**merged_data)
        # Filter the output values and return
        return {k: self._state[k] for k in self.outputs}
