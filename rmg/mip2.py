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

if args.verbose >2:
    if args.time:
        print("%12.3f " % (time() - all_start),end='')
    print("Load Instance data")
inst = Instance(args)

mip = mip.MIP2(instance=inst,args=args)
if args.verbose >0:
    if args.time:
        print("%12.3f " % (_time),end='')
    print("building MIP Model")

mip.build_model()

if args.dump:
    mip.write()

if args.verbose >0:
    if args.time:
        print("%12.3f " % (_time),end='')
    print("Solving Model")
solution = mip.solve()

print(' '.join(str(i) for i in inst.assign()))
print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
