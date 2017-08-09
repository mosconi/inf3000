#!/usr/bin/env python3

import argparse,os,sys
import numpy as np
from time import time

from rmg import roadef,common,instance,Instance,columngeneration

import rmg.columngeneration as CG

from gurobipy import tupledict

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser,
                                          columngeneration.parser
                                          
]
)


# Other parameters:

args = parser.parse_args()

if args.verbose is None:
    args.verbose = 1

if args.name:
    print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
    print(args)
    sys.exit(0)

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+17d" % x, })


inst = Instance(args)

cg = CG.CG4(instance=inst,args=args)

for ndeps in sorted(inst.sS):
    if args.verbose >3:
        print("building LP Model for services with %d deps" % ndeps)

    cg.build_lpmodel(ndeps)
    """
    for m in range(inst.nmach):
        if args.verbose >3:
            print(" building machine %5d (of %5d) MIP model " % (m,inst.nmach))
            cg.build_column_model(m,ndeps)
            if args.generate:
                if args.verbose >3:
                    print(" Pregenerate columns")
                pass
    """                
        
            

sys.exit(0)
if args.dump:
    cg.lpwrite(k=-1)

    
if args.verbose>3:
    print("Converting LP model to MIP model")
cg.lp2mip()

if args.dump:
    cg.write()
    
if args.verbose>3:
    print("Solve MIP model")

solution = cg.solve()

#print(solution)

if args.verbose>0:
    print(' '.join(str(i) for i in inst.assign()))
    print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
