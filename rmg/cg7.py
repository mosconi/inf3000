#!/usr/bin/env python3

import argparse,os,sys
import numpy as np
from time import time

from rmg import roadef,common,instance,columngeneration
from rmg.instance import Instance
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

if args.name:
    print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
    sys.exit(0)

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+20d" % x, })


inst = Instance(args)

cg = CG.CG3(instance=inst,args=args)

cg.build_lpmodel()

for m in range(inst.nmach):
    cg.build_column_model(m)

omega=0
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

alpha = .5

continue_cond = True

while continue_cond:
    cg.lpwrite()

    res = cg.solve_relax()
    
    if res is None:
        raise Exception("LP não ótimo")

    pi = alpha * best_pi + (1 - alpha) * res.pi
    mu = alpha * best_mu + (1 - alpha) * res.mu
    gamma = alpha * best_gamma + (1 - alpha) * res.gamma
    eta_lb = alpha * best_eta_lb + (1 - alpha) * res.eta_lb
    eta_ub = alpha * best_eta_ub + (1 - alpha) * res.eta_ub
    omikron_lb = alpha * best_omikron_lb + (1 - alpha) * res.omikron_lb
    omikron_ub = alpha * best_omikron_ub + (1 - alpha) * res.omikron_ub
    
    for m in range(inst.nmach):
        print(k,m)
        
        sigma = eta_lb[:,m] + eta_ub[:,inst.iN[m]] + omikron_lb[:,m] + omikron_ub[:,inst.iL[m]]

        cres = cg.compute_column(machine = m,
                                 pi = pi,
                                 gamma = gamma,
                                 sigma = sigma,
                                 mu = mu[m])

        if cres is None:
            raise Exception("problema ao calcular coluna")
        
        vres = cg.validate_column(cres,
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

        if ares.status != CG.CGAddStatus.Added:
            print(ares)
            raise Exception("problema ao adicionar coluna")


    if omega > - args.epslon:
        continue_cond = False
        break
        
    
    k+=1
    if k >= 4 and False:
        continue_cond = False
    
cg.lpwrite()

cg.lp2mip()

solution = cg.solve()

#print(solution)


print(' '.join(str(i) for i in solution.assign))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
