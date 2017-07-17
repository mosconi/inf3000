
from collections import namedtuple
from enum import Enum

RelaxSolution = namedtuple('RelaxSolution',['obj','pi','mu','eta_lb','eta_ub','gamma','omikron_lb','omikron_ub','rtime'])

CGColumn = namedtuple('CGColumn',['rc','obj','procs','g','servs','z','hsigma','ggamma','pixp','u','ut','a','d','b','pmc','mmc','rtime'])

CGValidate = namedtuple('CGValidate',['status'])

CGAdd = namedtuple('CGAdd',['status','var'])

CGSolution = namedtuple('CGSolution',['obj','X','assign'])

class CGValidateStatus(Enum):
    Valid = 1
    CalcMismatch = 2
    Invalid = 3

    

class CGAddStatus(Enum):
    Added = "A"
    Exist = "E"
    NotAdded = "N"
