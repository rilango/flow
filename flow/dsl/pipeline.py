import sys
import logging
import asyncio

from flow.config import Configuration
from flow.utils.stringutils import substitute
from .base import TaskBlock

log = logging.getLogger(__name__)


class Pipeline(TaskBlock):
    '''
    Used for defining pipelines.
    '''
    def __init__(self, **kwargs):
        TaskBlock.__init__(self, **kwargs)

    def __call__(self, *args, **kwargs):
        '''
        Parses the compute graph and retains it for later execution.
        '''
        assert len(args) > 0, 'Pipeline defined without tasks'

        log.debug(f'Processing pipeline "{self.name}"')
        TaskBlock.__call__(self, *args, **kwargs)

        self.inputs, self.outputs = self.effective_values()

        task_args = self.copy_attributes(**kwargs)
        task_args['pipeline_id'] = self.id

        return ExecutablePipeline(**task_args)


class ExecutablePipeline(Pipeline):

    async def __call__(self, **kwargs):
        inputs = {k: kwargs[k] for k in self.inputs if kwargs.get(k) is not None}
        assert len(self.inputs) == len(inputs), \
            f'Invalid # of inputs to "{self.name}" expected values for \"{", ".join(self.inputs)}\"'
        log.info(f'Running pipeline "{self.name}({self.id})"')
        log.debug(f'Pipeline arguments {kwargs}')

        # Update config with pipeline details.
        # TODO: Revisit for nested pipelines & loops. Hint platform.create_sub_ws
        config = Configuration()
        config.config['inputs'] = kwargs
        ws = config.platform.create_ws(uid=self.id, **config.config)
        config.config = {'pipeline': self.name,
                          'pipeline_id': self.id,
                          'pipeline_path': ws}

        log.info(f'Configuration {config.config}')
        self._state = config.platform.morph_arguments(config, **kwargs)

        log.info('Starting services...')
        await self._start_services()
        await self.execute_block()

    async def _start_services(self, **kwargs):
        # Start all services.
        for name, service in self._services.items():
            assert callable(service.func), f'Service "{name}" is not callable'
            from threading import Thread
            t = Thread(target = service)
            t.start()

        # Make sure all services are available.
        for name, service in self._services.items():
            retry_count = 0
            while retry_count < 30:
                if service.is_available():
                    log.info(f'Service found after {retry_count} retries.')
                    break
                else:
                    log.warning(f'Service not available. Retrying {retry_count}...')
                    import time
                    time.sleep(10)
                    retry_count += 1

    def _update_state(self, task, values):
        '''
        Updates the values to available inputs, based on task's output spec.
        '''
        if type(values) in [list, tuple, set]:
            # When multiple outputs are returned, merge to _state,
            # assuming the values are returned in the same order as in spec.
            assert len(task.outputs) == len(values), \
                f'Task "{task.name}" expected {len(task.outputs)} values, got {len(values)}'
            self._state.update(zip(task.outputs, values))

        elif type(values) in [dict]:
            # This is a special case for the loop condition. Just update after
            # checking the size of the dictionary.
            assert len(task.outputs) == len(values), \
                f'Task "{task.name}" expected {len(task.outputs)} values, got {len(values)}'
            self._state.update(values)

        else:
            # When single output is returned, we just update _state
            assert len(task.outputs) == 1, \
                f'Task "{task.name}" expected {len(task.outputs)} values, got 1'
            self._state[task.outputs[0]] = values

    async def _status_change(self, task, status, values):
        log.debug(f'Task {task.name} returned {values}')

        assert task.outputs is not None and values is not None, \
            f'Task "{task.name}" expected to return {", ".join(task.outputs)}, got None'

        log.debug(f'{task.name} returned {values}...')
        self._update_state(task, values)

        # At this point a task is completed and another could be waiting for
        # it to complete. TODO: Revisit: This is NOT a smart way to do this.
        if len(self._visited_tasks) < len(self._graph):
            await self.execute_block()

        if status == 'completed':
            self._completed_tasks.add(task.name)

    async def _execute_task(self, task, inputs):
        inputs = inputs.copy()
        inputs.update(self._state)
        if task.constants is not None:
            inputs.update(task.constants)

        log.debug(f'Calling "{task.name}" with inputs {inputs}...')
        result = await task(**inputs)
        await self._status_change(task, 'completed', result)

    def _construct_input(self, task):
        '''
        Constructs the input dictionary for a task. Inputspec is used
        to get values from _state.
        '''
        if task.inputs is None or len(task.inputs) == 0:
            return {}
        try:
            inputs = {d: self._state[d] for d in task.inputs}
            if len(inputs) == len(task.inputs):
                return inputs
            else:
                log.debug(f'Insufficent inputs for task "{task.name}"')
                return None
        except KeyError:
            return None

    async def execute_block(self):
        '''
        Executes all tasks that fulfill all pre-requisites.

        Returns:
            True if a task is executed, False if none met the pre-req.
        '''

        # Identify tasks for which all inputs are available
        tasks = []
        for task in self._graph.values():
            assert callable(task), f'Task "{task.name}" is not callable'
            if task.name in self._visited_tasks:
                continue

            inputs = self._construct_input(task)
            if inputs is not None:
                # TODO: Revisit for throttling and resource management
                log.debug(f'Executing task {task.name} with inputs {inputs}...')
                self._visited_tasks.add(task.name)

                tasks.append(asyncio.create_task(
                    self._execute_task(task, inputs)))

        await asyncio.gather(*tasks)
        return len(tasks) > 0
