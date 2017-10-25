#!/usr/bin/env python

import argparse,os,sys
import numpy as np

np.set_printoptions(formatter={'float_kind': lambda x: "%+17.3f" %
x,'int': lambda x: "%+17d" % x, })

from time import time

from rmg import roadef,common,instance,Instance,ValidateStatus
import rmg.SAT as SAT

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser
]
)

args = parser.parse_args()

if args.verbose is None:
    args.verbose = 1

all_start =0

def msg(level, msg, timeref=None):
    global all_start 

    if args.verbose >level:
        if args.time:
            _time = time()
            print("%12.3f " % (_time  - all_start),end='')
            if timeref is not None:
                print("%12.3f " % (_time - timeref),end='')
                
        print(msg)


def main():
    
    if args.name:
        print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
        return


    
    global all_start 
    all_start = time()
    msg(0,"Load Instance data")
    #msg(2,"Load Instance data")
            
    inst = Instance(args)
    msg(0,"done")

    msg(0,"sat class")
    sat = SAT.SAT1(inst)
    msg(0,"build model")
    sat.build_model()
    msg(0,"done")

    msg(0,"Solve until UNSAT")

    
    best_assign = None
    best_obj = np.iinfo(np.int32).max

    k= 0 
    while True:

        solution = sat.solve()
        
        if not solution.status:
            break

        invalid_column = False
        for m in inst.M:
            v = inst.mach_validate(m,solution.assignment[:,m])
            if v.status != ValidateStatus.Valid:
                msg(1,"invalid col for mach %d @ iter %d %s" %(m,k,v.status))
                sat.exclude_column(m,solution.assignment[:,m])
                invalid_column = True
        k+=1
        if invalid_column:
            continue

        v = inst.map_validate(solution.assignment)

        msg(2,"assign is %s" %(v.status))
        if v.status == ValidateStatus.Valid and best_obj > v.obj:
            msg(2,"improving obj %d -> %d" %(best_obj, v.obj))
            best_obj = v.obj
            best_assign = solution.assignment

        sat.exclude_solution(solution)
        
    msg(0,"done")
    
    msg(0,"sol")

    print(best_obj)
    print(best_assign)
    msg(0,"done")
