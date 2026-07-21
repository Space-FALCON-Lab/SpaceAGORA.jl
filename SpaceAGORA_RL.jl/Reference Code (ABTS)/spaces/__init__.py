from spaces.space import Space
from spaces.box import Box
from spaces.discrete import Discrete
from spaces.multi_discrete import MultiDiscrete
from spaces.multi_binary import MultiBinary
from spaces.tuple import Tuple
from spaces.dict import Dict

from spaces.utils import flatdim
from spaces.utils import flatten_space
from spaces.utils import flatten
from spaces.utils import unflatten

__all__ = ["Space", "Box", "Discrete", "MultiDiscrete", "MultiBinary", "Tuple", "Dict", "flatdim", "flatten_space", "flatten", "unflatten"]