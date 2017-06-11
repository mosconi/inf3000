#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from instance import Instance

from gurobipy import *

from columngeneration import CG5

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+20.3f" % x,'int': lambda x: "%+20d" % x, })

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=True
epslon = 10**-1
xi = 10**3
hsz = 10
slack = 0.02

try:
    opts, args = getopt.getopt(sys.argv[1:],"tqsn:m:a:o:h",["tex","quiet","save","name=","model=","assign=","output=","help"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit()
    
for o,a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ( "-a", "--assign"):
        assignfile = a
    elif o in ( "-m", "--model"):
        modelfile = a
    elif o in ( "-n", "--name"):
        name = a
    elif o in ( "-o", "--output"):
        outputfile = a
    elif o in ( "-s", "--save"):
        savemodel=True
    elif o in ( "-q", "--quiet"):
        verbose=False
    elif o in ( "-t", "--tex"):
        tex=True
    else:
        assert False, "unhandled option"
    
if modelfile is None:
    sys.exit("Model file not defined")

if assignfile is None:
    sys.exit("Assign file not defined")

if not os.path.exists(modelfile):
    sys.exit("File %s not found" % modelfile)

if not os.path.exists(assignfile):
    sys.exit("File %s not found" % assignfile)


inst = Instance(model=modelfile, assign=assignfile)

cg = CG5(inst, epslon=0.5)


# Initialize restricted master problem
cg.build_model()

# Initialize vars

k=0
impr=0

continue_cond = True

best_pi = np.zeros(inst.nproc)
best_alpha = np.zeros(inst.nmach)
best_omega = - np.inf

beta0 = .5
betaJ = .95

beta = beta0

from multiprocessing import Pool
from functools import partial
pool = Pool(processes=4)

def comp_m(mdl,k,m,pi,alpha):
    print(" iter %05d  mach %05d => " % (k,m),end='')
    model = mdl.copy()
    (roadef, w , q,model)= model.compute_column(m, pi, alpha[m])
    
    print("%20.3f (roadef: %20.3f rat: %20.10f)" % (w, roadef,abs(w/(roadef+1))))
    model.validate_column(q,m,w,roadef,pi,alpha[m],epslon,model)
    if len(q) < 100:
        print("add col", q)
        
    mdl.mip_add_col(obj=roadef, col =q, machine=m)
        
    #cols = cg.generate_companion_columns(machine=m,col=q,obj=roadef
    #for (_m,_obj,_col) in cols
    return w

while continue_cond:
    # solve RDWM 
    (z_rm, pi_rm, alpha_rm) = cg.solve_relax(k=k)

    if np.isnan(z_rm):
        break

    # Compute pi_stabilized
    pi = beta * best_pi + ( 1 - beta ) * pi_rm
    alpha = beta * best_alpha  + ( 1 - beta ) * alpha_rm

    # omega é a soma dos custos reduzidos dos sub problemas
    omega = 0
    print(" iter %05d  obj: %+20.3f" %(k, z_rm) )

    # solve L(pi_st)
#    for m in range(inst.nmach):
#        print(" iter %05d  mach %05d => " % (k,m),end='')
#        (roadef, w , q,model)= cg.compute_column(m, pi, alpha[m])
##        if w > - epslon:
##            print("for break")
##            break

#        print("%20.3f (roadef: %20.3f rat: %20.10f)" % (w, roadef,abs(w/(roadef+1))))
#        cg.validate_column(q,m,w,roadef,pi,alpha[m],epslon,model)
#        omega += w
#        if len(q) < 100:
#            print("add col", q)
#        cg.mip_add_col(obj=roadef, col =q, machine=m)

        #cols = cg.generate_companion_columns(machine=m,col=q,obj=roadef
        #for (_m,_obj,_col) in cols

    _w = pool.map(partial(comp_m,mdl=cg,k=k,pi=pi,alpha=alpha),[m for m in range(inst.nmach)])

    print(_w)

    print("omega: %20.3f (%20.3f) sum: %20.3f (%20.3f) N: %d, beta: %0.3f)" % (omega, best_omega, z_rm + omega,z_rm + best_omega, impr, beta))

    
    if omega > best_omega:
        impr += 1
        beta = max(.5, beta0*(betaJ**int(impr/10)))

        print("best omega improvement %20.3f -> %20.3f (%d, %0.3f)" % (best_omega, omega, impr, beta))
        best_omega = omega
        best_pi = pi
        best_alpha = alpha
        #input("Press Enter to continue...")

    if  best_omega >- epslon:
        contine_cond = False
        break
        
    k+=1
#    if k>= 250:
#        break

(obj, X, alloc) = cg.solve_mip()

print(X)
