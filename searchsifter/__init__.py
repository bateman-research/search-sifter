from .Family import Family
from . import relationships
from . import scripts

import logging
import sys

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
