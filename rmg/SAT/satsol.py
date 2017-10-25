
from collections import namedtuple
from enum import Enum

SATSolution = namedtuple('SATSolution',['status','assignment'])

class SATStatus(Enum):
    Valid = 1
    Invalid = 2
