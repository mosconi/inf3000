#!/usr/bin/env python3

import argparse,os,sys
import numpy as np
from time import time

from rmg import roadef,common,instance,mip
from rmg.instance import Instance

from gurobipy import *

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser,
                                          mip.parser
                                          
]
)


# Other parameters:

args = parser.parse_args()
if args.name:
    print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
    sys.exit(0)

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+20d" % x, })


inst = Instance(args)

mip = mip.MIP2(instance=inst,args=args)

mip.build_model()
mip.write()

solution = mip.solve()

print(' '.join(str(i) for i in inst.assign()))
print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
