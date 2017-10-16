#!/usr/bin/env python3

import argparse,os,sys
import numpy as np

np.set_printoptions(formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+17d" % x, })

from time import time

from rmg import roadef,common,instance,Instance

from satispy import Variable,Cnf
from satispy.io import DimacsCnf
from satispy.solver import Minisat


from itertools import combinations

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

def sat(instance):
    expr = Cnf()

    nproc = instance.nproc
    nmach = instance.nmach
    P = instance.P
    M = instance.M
    S = instance.S
    N = instance.N
    L = instance.L

    print("creating vars - x")
    x = np.array([[Variable('x[%d,%d]' % (p,m)) for m in M] for p in P])
    print("creating vars - h")
    h = np.array([[Variable('h[%d,%d]' % (s,n)) for n in N] for s in S])
    print("creating vars - o")
    o = np.array([[Variable('o[%d,%d]' % (s,l)) for l in L] for s in S])

    print("H[s,n]")
    for s in S:
        for n in N:
            pres_expr = Cnf()
            for p in S[s]:
                for m in N[n]:
                    pres_expr |=x[p,m]
                    
            expr &= h[s,n] >> pres_expr
            expr &= pres_expr >> h[s,n]
            for sd in instance.sdep[s]:
                expr &= h[s,n] >> h[sd,n]
    print("O[s,l]")
    for s in S:
        for l in L:
            pres_expr = Cnf()
            for p in S[s]:
                for m in L[l]:
                    pres_expr |=x[p,m]
                    
            expr &= o[s,l] >> pres_expr
            expr &= pres_expr >> o[s,l]
    
    print("X[p,m]")
    
    for p in P:
        p_constr1 = Cnf()
        p_constr2 = Cnf()
        for m in M:
            p_constr1 |= x[p,m]
        for _m in combinations(M,2):
            p_constr2 &= x[p,_m[0]] >> -x[p,_m[1]]
        expr &= p_constr1 & p_constr2

    print("X[p,m] - conflito")
        
    for m in M:
        for s in S:
            conf_constr = Cnf()
            if len(S[s]) ==1: continue
            for _p in combinations(S[s],2):
                conf_constr &= x[_p[0],m] >> -x[_p[1],m]
            expr &= conf_constr


    io = DimacsCnf()
    s = io.tostring(expr)
    print(s)
    solver = Minisat('minisat %s %s')

    solution = solver.solve(expr)

    if not solution.success:
        print("sem solucao")
        exit(0)

    xsol = np.empty((nproc,nmach),dtype=np.int32)
    for p in P:
        for m in M:
            print(p,m,solution[x[p,m]])
            print(io.varname(x[p,m]))
            xsol[p,m] = solution[x[p,m]]

    print(xsol)
    print(np.all(xsol.sum(axis=1)==1))
    for m in M:
        print(instance.mach_validate(m,xsol[:,m]))
    
    
    
    

def main():
    
    if args.name:
        print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
        return

    global all_start
    all_start = time()
    print("Load Instance data")
    #msg(2,"Load Instance data")
            
    inst = Instance(args)
    print("done")

    sat(inst)



