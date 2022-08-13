from distutils.util import get_platform
from multiprocessing.sharedctypes import Value
from .container import ContainerizedTask, ContainerPlatform



class PlatformBuilder():

    @staticmethod
    def get_platform(**kwargs):
        '''
        Creates the desired platform implementation.

        :param str platform_type: Type of platform. Currently supports 'containers'
        '''
        if kwargs['platform_type'] == 'containers':
            return ContainerPlatform(**kwargs)
        else:
            raise ValueError(
                f'Unsupported platform "{kwargs["platform_type"]}"')
