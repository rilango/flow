from subprocess import call
import uuid
import logging

log = logging.getLogger(__name__)


class Component(object):

    def __init__(self, **kwargs):
        # self.id is a unique identifier and is tied to the implementation of
        # the task. For example, Platform plugins will use this to create a
        # unique location/workspace.
        self.id = kwargs.get('id')
        if self.id is None:
            self.id = str(uuid.uuid1())

        self.name = kwargs.get('name')
        self.inputs = kwargs.get('inputs')
        self.outputs = kwargs.get('outputs')
        self.constants = kwargs.get('constants')
        self.params = kwargs.get('params')
        self._service = kwargs.get('_service')

        # This is a list of string representing namespace of the task. Each
        # instance of task has to be represented uisng a unique path. Which
        # needs additional care during runtime.
        self._node_path = [self.id]

        assert self.name is not None, 'Name is required'
        if self.inputs is not None and type(self.inputs) is not list:
            self.inputs = [self.inputs]

        if self.constants is not None and type(self.constants) is not dict:
            raise TypeError(f'"constants" in {self.name} must be a dictionary')

    @property
    def node_path(self):
        return self._node_path

    @node_path.setter
    def node_path(self, parent_uid):
        if parent_uid is None:
            return
        self._node_path.insert(0, parent_uid)

        # Update the node path of all the children
        if hasattr(self, '_graph'):
            for task in self._graph.values():
                task.node_path = parent_uid

    def copy_attributes(self, **kwargs):
        task_args = {}
        task_args.update(vars(self))
        task_args.update(kwargs)

        return task_args


class BaseTask(Component):
    '''
    Base abstract class for implementing tasks. It defines the structure to be
    used for defining a task.

    :param name: Name of the task.
    :param inputs: List of inputs to the task.
    :param outputs: List of outputs from the task.

    '''

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for a task.
        '''
        Component.__init__(self, **kwargs)

    def __call__(self, **kwargs):
        pass

class TaskBlock(Component):
    '''
    Abstract class for implementing composite structures consisting of one or
    more tasks.
    '''

    def __init__(self, **kwargs):
        '''
        Reads the minimum required parameters for a loop.
        '''
        Component.__init__(self, **kwargs)
        # Output is the combination of all outputs of the tasks
        self.outputs = None

        # Tasks are managed at TaskBlock level to ensure scope of all variable
        # (input/output) are maintained at that level. Following variables are
        # used to maintain the state of the TaskBlock.
        self._graph = kwargs.get('_graph', {})
        self._state = kwargs.get('_state', {})
        self._visited_tasks = kwargs.get('_visited_tasks', set())
        self._completed_tasks = kwargs.get('_completed_tasks', set())

        # To maintain a list of all services. This will be used for cleanup
        self._services = kwargs.get('_services', {})

    def __call__(self, *args, **kwargs):
        '''
        Parses the compute graph and retains it for later execution.
        '''
        assert len(args) > 0, 'Pipeline defined without tasks'

        for task in args:
            assert callable(task), f'{task.name} is not callable'
            if task.name in self._graph:
                raise ValueError(f'Task "{task.name}" already defined')

            if task._service:
                self._services[task.name] = task
            else:
                task.node_path = self.id
                self._graph[task.name] = task

    def effective_values(self):
        '''
        Identifies required inputs required for a block to start and the output
        values generated within the block.
        '''
        # Compute effective inputs
        all_inputs = set()
        all_outputs = set()

        for task in self._graph.values():
            log.debug(f'{task.name}(task.inputs)->{task.outputs}')

            if task.inputs is not None:
                # Ignore internal inputs. They start with _. i.e. _loop_count
                all_inputs.update([k for k in task.inputs if not k.startswith('_')])

            if task.outputs is not None:
                outputs = set(task.outputs) - set(task.inputs)
                all_outputs.update(outputs)

        all_inputs = list(all_inputs - all_outputs)

        assert len(all_inputs) > 0, 'No inputs found'
        return all_inputs, list(all_outputs)
