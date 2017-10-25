#!/usr/bin/env python3

import argparse,os,sys
import numpy as np

np.set_printoptions(formatter={'float_kind': lambda x: "%+17.3f" %
x,'int': lambda x: "%+17d" % x, })

from time import time

from rmg import roadef,common,instance,Instance

from satispy import Variable,Cnf
from satispy.io import DimacsCnf
from satispy.solver import Minisat

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
    nserv = instance.nserv
    nneigh = len(instance.N)
    nloc = len(instance.L)

    P = instance.P
    M = instance.M
    S = instance.S
    N = instance.N
    L = instance.L

    print(" %10.3f creating vars - x" % (time() - all_start))
    x = np.empty((nproc,nmach), dtype=object)
    for p in P:
        for m in M:
            x[p,m]=Variable('x[%d,%d]' % (p,m))
            
    print(" %10.3f creating vars - h" % (time() - all_start))
    h = np.empty((nserv,nneigh), dtype=object)
    for s in S:
        for n in N:
            h[s,n] = Variable('h[%d,%d]' % (s,n))

    print(" %10.3f creating vars - o" %(time() - all_start))
    o = np.empty((nserv,nloc), dtype=object)
    for s in S:
        for l in L:
            o[s,l] = Variable('o[%d,%d]' % (s,l))

    print(" %10.3f H[s,n]" %( time() - all_start))
    for s in S:
        #print(" %10.3f H[%6d,n]" %( time() - all_start,s))
        for n in N:
            pres_expr = Cnf()
            for p in S[s]:
                for m in N[n]:
                    pres_expr |= x[p,m]
                    
            expr &= h[s,n] >> pres_expr
            expr &= pres_expr >> h[s,n]
            for sd in instance.sdep[s]:
                expr &= h[s,n] >> h[sd,n]
                
    print(" %10.3f O[s,l]" %( time() - all_start))
    for s in S:
        #print(" %10.3f O[%6d,l]" %( time() - all_start,s))
        if instance.delta[s] ==1 : continue
        for l in L:
            pres_expr = Cnf()
            for p in S[s]:
                for m in L[l]:
                    pres_expr |=x[p,m]
                    
            expr &= o[s,l] >> pres_expr
            expr &= pres_expr >> o[s,l]
   
    print(" %10.3f X[p,m]" %( time() - all_start))
    
    for p in P:
        #print(" %10.3f X[%6d, m]" %( time() - all_start, p))
        p_constr1 = Cnf()
        p_constr2 = Cnf()
        for m in M:
            p_constr1 |= x[p,m]

            for m2 in range(m+1,nmach):
                p_constr2 &= x[p,m] >> -x[p,m2]
        expr &= p_constr1 & p_constr2

    print(" %10.3f X[p,m] - conflito" %(time() - all_start))
        
    for m in M:
        for s in S:
            conf_constr = Cnf()
            if len(S[s]) == 1: continue
            for i1 in range(len(S[s])):
                for i2 in range(i1+1,len(S[s])):
                    conf_constr &= x[S[s][i1],m] >> -x[S[s][i2],m]
            expr &= conf_constr


    print(" %10.3f solving" %( time() - all_start))

    io = DimacsCnf()
    #s = io.tostring(expr)
    #print(s)
    solver = Minisat('minisat %s %s')

    solution = solver.solve(expr)
    print(" %10.3f done" %( time() - all_start))

    if not solution.success:
        print(" %10.3f sem solucao" % (time() - all_start))
        exit(0)

    print(" %10.3f solucao" % (time() - all_start))

    xsol = np.empty((nproc,nmach),dtype=np.int32)
    for p in P:
        for m in M:
            #print(p,m,solution[x[p,m]])
            #print(io.varname(x[p,m]))
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
    print(" %10.3f Load Instance data" % (time() - all_start))
    #msg(2,"Load Instance data")
            
    inst = Instance(args)
    print(" %10.3f done" % (time() - all_start))

    sat(inst)
    print(" %10.3f done" % (time() - all_start))



