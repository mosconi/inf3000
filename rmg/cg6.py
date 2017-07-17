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

cg.build_lpmodel()

for m in range(inst.nmach):
    cg.build_column_model(m)

k=0
impr=0
best_zrm = + np.inf
best_omega = - np.inf
best_pi = np.zeros(inst.nproc,dtype=np.float64)                        # Duais de processos
best_mu = np.zeros(inst.nmach,dtype=np.float64)                        # Duais de máquina
best_gamma = np.zeros(inst.nserv,dtype=np.float64)                     # Duais de serviço migrado
best_eta_lb = np.zeros((inst.nserv,inst.nmach),dtype=np.float64)       # Dual de h[s,n,m]
best_eta_ub = np.zeros((inst.nserv,len(inst.N)),dtype=np.float64)      # Dual de h[s,n]
best_omikron_lb = np.zeros((inst.nserv,inst.nmach),dtype=np.float64)   # Dual de o[s,l,m]
best_omikron_ub = np.zeros((inst.nserv,len(inst.L)),dtype=np.float64)  # Dual de o[s,l]

alpha = 0

continue_cond = True

res = cg.solve_relax()
first_obj = res.obj


while continue_cond:
    if args.verbose >1:
        print("-"*78)
        print("%5d master" % (k),end=' ',flush=True)

    r_start = time()
    if args.dump:
        cg.lpwrite(k)

    res = cg.solve_relax()
    r_stop = time()

    if args.verbose >1:
        print("%20.3f %20.3f %8.3f (%8.3f) %20.3f" % (res.obj,first_obj, r_stop - r_start, res.rtime, res.obj - first_obj))
    
    if res is None:
        raise Exception("LP não ótimo")

    pi = alpha * best_pi + (1 - alpha) * res.pi
    mu = alpha * best_mu + (1 - alpha) * res.mu
    gamma = alpha * best_gamma + (1 - alpha) * res.gamma
    eta_lb = alpha * best_eta_lb + (1 - alpha) * res.eta_lb
    eta_ub = alpha * best_eta_ub + (1 - alpha) * res.eta_ub
    omikron_lb = alpha * best_omikron_lb + (1 - alpha) * res.omikron_lb
    omikron_ub = alpha * best_omikron_ub + (1 - alpha) * res.omikron_ub

    omega = 0
    add_col = False

    for m in range(inst.nmach):
        if args.verbose>1:
            print("%5d  %5d" % (k,m),end=' ',flush=True)

        c_start = time()
        sigma = eta_lb[:,m] + eta_ub[:,inst.iN[m]] + omikron_lb[:,m] + omikron_ub[:,inst.iL[m]]

        cg.column_prepare(machine = m,
                          pi = pi,
                          gamma = gamma,
                          sigma = sigma,
                          mu = mu[m])
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
                                      pi = pi,
                                      gamma = gamma,
                                      sigma = sigma,
                                      mu = mu[m])

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
        print("%5d  omega %20.3f %20.3f %8.3f (%8.3f)" %(k,omega,best_omega, time() - r_start , 0),end=' ')
        
    if omega > - args.epslon or not add_col:
        continue_cond = False

    k+=1
    if k >= 4 and False:
        continue_cond = False

    if continue_cond:
        print("C")
    else:
        print("T")
        
cg.lp2mip()

solution = cg.solve()

#print(solution)

if args.verbose>0:
    print(' '.join(str(i) for i in inst.assign()))
    print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
