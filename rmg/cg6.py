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

cg = CG.CG2(instance=inst,args=args)

if args.verbose >3:
    print("building LP Model")

cg.build_lpmodel()

for m in range(inst.nmach):
    if args.verbose >3:
        print(" building machine %5d (of %5d) MIP model" % (m,inst.nmach))
    cg.build_column_model(m)
    if args.generate:
        for p in range(inst.nproc):
            m_assign = inst.mach_map_assign(m).copy()
            m_assign[p] = not m_assign[p]
            if inst.mach_validate(machine=m,map_assign=m_assign):
                pass


            
continue_cond = True

res = cg.solve_relax()
first_obj = res.obj

stab = CG.Stabilization(inst, args)

alpha = stab.alpha0()

k = 0
while continue_cond:
    if args.verbose >1:
        print("-"*(int(columns)-2))
        print("%5d master" % (k),end=' ',flush=True)

    r_start = time()
    if args.dump:
        cg.lpwrite(k)

    res = cg.solve_relax()
    if args.dump:
        cg.lpwritesol(k)

    r_stop = time()

    if args.verbose >1:
        print("%20.3f %20.3f %8.3f (%8.3f) %20.3f" % (res.obj,first_obj, r_stop - r_start, res.rtime, res.obj - first_obj))
    
    if res is None:
        raise Exception("LP não ótimo")
    

    stduals = stab.alpha(res)
    omega = 0
    add_col = False

    for m in range(inst.nmach):
        if args.verbose>1:
            print("%5d  %5d" % (k,m),end=' ',flush=True)

        c_start = time()
        sigma = stduals.eta_lb[:,m] + stduals.eta_ub[:,inst.iN[m]] + stduals.omikron_lb[:,m] + stduals.omikron_ub[:,inst.iL[m]]

        cg.column_prepare(machine = m,
                          pi = stduals.pi,
                          gamma = stduals.gamma,
                          sigma = sigma,
                          mu = stduals.mu[m])
        if args.dump:
            cg.column_write(machine = m, k = k)
        
        cres = cg.column_compute(machine = m)

        c_stop = time()
        if cres is None:
            raise Exception("problema ao calcular coluna")


        if args.dump:
            cg.column_writesol(machine = m, k = k)

        if args.validate:
            vres = cg.column_validate(cres,
                                      machine = m,
                                      pi = stduals.pi,
                                      gamma = stduals.gamma,
                                      sigma = sigma,
                                      mu = stduals.mu[m])

            if vres.status != CG.CGValidateStatus.Valid:
                print(vres)
                raise Exception("problema ao validar coluna")
        
        omega += cres.rc
        ares = cg.lp_add_col(m,cres)

        if ares.status == CG.CGAddStatus.Added:
            add_col = True
            
        elif ares.status == CG.CGAddStatus.NotAdded:
            print(ares)
            raise Exception("problema ao adicionar coluna")
        
        if args.verbose>1:
            print("%20.3f %20.3f %8.3f (%8.3f) %s" %(cres.rc, cres.obj, c_stop -c_start ,cres.rtime ,ares.status.value ) )
        
    if args.verbose>1:
        print("%5d  omega %20.3f %20.3f %8.3f (%8.3f)" %(k,omega,stab.best_omega(), time() - r_start , 0),end=' ')

    
    if omega > - args.epslon or not add_col:
        continue_cond = False

    k+=1
    if k >= 4 and False:
        continue_cond = False

    if continue_cond:
        print("C %0.6f %d" %(alpha, stab.improvements()))
    else:
        print("T")

    alpha = stab.compute(omega = omega,stabdual = stduals)


if args.verbose>3:
    print("Converting LP model to MIP model")
cg.lp2mip()

if args.verbose>3:
    print("Solve MIP model")

solution = cg.solve()

#print(solution)

if args.verbose>0:
    print(' '.join(str(i) for i in inst.assign()))
    print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
