#!/usr/bin/env python3

import argparse,os,sys
import numpy as np
from time import time

import roadef
import common
import instance
import mip 

from gurobipy import *

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser,
                                          mip.parser]
)


# Other parameters:

args = parser.parse_args()

print(args)

print(roadef.args)

if roadef.args.name:
    print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
    sys.exit(0)

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+20d" % x, })


inst = instance.Instance(opened=True,model=args.instance_filename,
                assign=args.original_solution_filename)

mip = mip.MIP1(instance=inst)

mip.build_model()

(obj,X,solution) = mip.solve()


#if args.new_solution_filename:
#    args.new_solution_filename.write(' '.join(str(i) for i in solution))
    
