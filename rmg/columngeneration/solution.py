
from collections import namedtuple

RelaxSolution = namedtuple('RelaxSolution',['obj','pi','mu','eta_lb','eta_ub','gamma','omikron_lb','omikron_ub'])

CGColumn = namedtuple('CGColumn',['rc','obj','procs','g','servs'])
