import io
import pytest
import logging
import asyncio

from flow.dsl import task, pipeline, loop, condition, service, when
from flow.config import Configuration

log = logging.getLogger(__name__)


def test_create_task():
    # Create a valid Task
    task(name='test', inputs='s', outputs='value')(print)

    # Create a Task without outputs
    with pytest.raises(AssertionError, match=r".*outputs.*"):
        task(name='test', inputs='s')(print)

    # Create a Task without input
    with pytest.raises(AssertionError, match=r".*must have inputs.*"):
        task(name='test')(print)


def test_create_service():
    # Create a valid Service
    service(name='test')(print)

    # Create a Service with inputs
    with pytest.raises(AssertionError, match=r".*must not have.*"):
        service(name='test', inputs='s')(print)

    # Create a Service with outputs
    with pytest.raises(AssertionError, match=r".*must not have.*"):
        service(name='test', outputs='s')(print)


def test_create_when():
    # Create a loop terminal condition
    when(name='check score')(
        task(name='test1', inputs='s', outputs='value')(print),
        task(name='test2', inputs='s', outputs='value')(print),
        condition(name='check score', inputs='s', outputs='value')(print),
    )

    # Create a lop without terminal condition
    with pytest.raises(AssertionError, match=r".*must contain.*"):
        when(name='check score')(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test2', inputs='s', outputs='value')(print)
        )

    # Create a loop with duplicate name
    with pytest.raises(ValueError, match=r".*already.*"):
        when(name='Loop1')(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test1', inputs='s', outputs='value')(print),
            condition(name='Loop1', inputs='s', outputs='value')(print),
        )

    # Empty loop
    with pytest.raises(AssertionError, match=r".*at least one.*"):
        when(name='Empty Loop')(
            condition(name='Empty Loop', inputs='s', outputs='value')(print),
        )

def test_create_loop():
    # Create a loop terminal condition
    loop(name='check score')(
        task(name='test1', inputs='s', outputs='value')(print),
        task(name='test2', inputs='s', outputs='value')(print),
        condition(name='check score', inputs='s', outputs='value')(print),
    )

    # Create a lop without terminal condition
    with pytest.raises(AssertionError, match=r".*must contain.*"):
        loop(name='check score')(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test2', inputs='s', outputs='value')(print)
        )

    # Create a loop with duplicate name
    with pytest.raises(ValueError, match=r".*already.*"):
        loop(name='Loop1')(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test1', inputs='s', outputs='value')(print),
            condition(name='Loop1', inputs='s', outputs='value')(print),
        )

    # Empty loop
    with pytest.raises(AssertionError, match=r".*at least one.*"):
        loop(name='Empty Loop')(
            condition(name='Empty Loop', inputs='s', outputs='value')(print),
        )


def test_create_pipeline():
    # Create a valid Pipeline
    pipeline(name='virtual_screening')(
        task(name='test1', inputs='s', outputs='value')(print),
        task(name='test2', inputs='s', outputs='value')(print)
    )

    # Create a Pipeline without name
    with pytest.raises(AssertionError, match=r".*required.*"):
        pipeline()(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test2', inputs='s', outputs='value')(print)
        )

    # Create a Pipeline without tasks
    with pytest.raises(AssertionError, match=r".*tasks.*"):
        pipeline(name='virtual_screening')()

    # Create tasks with duplicate names
    with pytest.raises(ValueError, match=r".*already.*"):
        pipeline(name='virtual_screening')(
            task(name='test1', inputs='s', outputs='value')(print),
            task(name='test1', inputs='s', outputs='value')(print)
        )


def returnOne(**kwargs):
    return 1


def returnTwo(**kwargs):
    return 1, 1

def generation(**kwargs):
    kwargs['generation'] = kwargs['generation'] + 1
    return kwargs

def condition_for_loop(**kwargs):
    log.info(f'Executing returnTrue input {kwargs}')
    return kwargs['_loop_count'] < 5

def condition_for_when(**kwargs):
    log.info(f'Executing returnTrue input {kwargs}')
    return True

def create_config():
    config = '''
    platform_type: containers
    platform_ws_root: /tmp/shared
    platform_db: /data/db/embedding_cache_cddd.sqlite3

    '''
    data = io.StringIO()
    data.write(config)
    config = Configuration()
    config.load_config(data)
    log.info(f'Starting pipeline with config: {config.config}')
    return config

def test_pipeline():

    config = create_config()
    # Create a valid Pipeline, This pipeline will expect a, b and d as inputs
    vs_pipline =\
        pipeline(name='virtual_screening')(
            task(name='ProteinPreparation',
                 inputs=['a', 'b'],
                 outputs='c')(returnOne),
            task(name='GenerateMolecule',
                 inputs=['b', 'c', 'd'],
                 outputs='generation')(returnOne),

            loop(name='Check Score', increment=generation)(
                task(name='GenerateConformers',
                     inputs=['generation'],
                     outputs='f')(returnOne),
                task(name='PrepareLigands',
                     inputs='f',
                     outputs='g')(returnOne),
                task(name='GeneratePose',
                     inputs='g',
                     outputs='h')(returnOne),
                when(name='GenerateWhen')(
                    condition(name='GenerateWhen',
                              inputs='generation')(condition_for_when),
                    task(name='GenerateMolecule1',
                         inputs=['h'],
                         outputs=['z'])(returnOne),
                ),
                condition(name='Check Score',
                          inputs='_loop_count',
                          outputs='z')(condition_for_loop),
            )
        )

    log.info(vs_pipline)
    assert len(vs_pipline._graph) == 3, 'Two tasks and one loops i.e. 3 nodes'
    assert 'Check Score' in vs_pipline._graph, 'check score loop not found'

    assert len(vs_pipline._graph['Check Score']._graph) == 4, \
        'Four tasks and the condition should get absorbed by the loop'
    assert vs_pipline._graph['Check Score'].condition.name == vs_pipline._graph['Check Score'].name, \
        'Loop did not absorb the right condition'

    assert len(vs_pipline.inputs) == 3, \
        f'Expected 4({vs_pipline.inputs}) inputs got {len(vs_pipline.inputs)}({vs_pipline.inputs})'

    log.info(f'Pipeline {vs_pipline.name} needs {vs_pipline.inputs}')

    inputs = {'a': 10, 'b': 20, 'd': 30}
    asyncio.run(vs_pipline(**inputs))

    log.info(config.config)


def test_pipeline_too_many_returns():
    config = create_config()
    # Create a valid Pipeline, This pipeline will expect a, b and d as inputs
    vs_pipline = pipeline(name='virtual_screening')(
        task(name='ProteinPreparation', inputs=[
             'a', 'b'], outputs='c')(returnTwo),
    )

    inputs = {'a': 10, 'b': 20}
    with pytest.raises(AssertionError, match=r".*expected.*"):
        asyncio.run(vs_pipline(**inputs))
