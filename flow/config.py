import yaml
import logging
from functools import singledispatch
from io import StringIO

from flow.utils.singleton import Singleton
from flow.platform import PlatformBuilder


log = logging.getLogger(__name__)


@singledispatch
def read_config(config):
    return NotImplemented


@read_config.register(StringIO)
def _(config):
    return yaml.safe_load(config.getvalue())


@read_config.register(str)
def _(config):
    with open(config, 'r') as stream:
        return yaml.safe_load(stream)


class Configuration(object, metaclass=Singleton):
    '''
    Encapsulates environment related properties that a pipeline may need.
    '''

    def __init__(self):
        self._config = {}
        self.platform = None

    def load_config(self, config_file:str):
        '''
        Load configuration. The config should contain minimul information
        required by the target platform.
        '''
        pipeline_config = read_config(config_file)
        self.update(**pipeline_config)

        self.platform = PlatformBuilder.get_platform(**self._config)

    def get(self, property_name, default=None):
        """
        Returns values from local configuration.
        """
        value = self._config.get(property_name)
        if value is None:
            log.warn(f'{property_name} not found, returing default.')
            return default
        return value

    def update(self, **kwargs):
        self.config = kwargs

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, values: dict):
        '''
        This setter does not overwrite the available inputs. It only updates.
        '''
        if values is not None:
            self._config.update(values)

    def get_pipeline_id(self):
        return self.get('pipeline_id')

    def get_pipeline_path(self):
        return self.get('pipeline_path')
