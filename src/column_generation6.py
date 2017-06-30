#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from instance import Instance

from gurobipy import *

from columngeneration import CG5

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+20d" % x, })

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name=None
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

print("Load instance")
inst = Instance(model=modelfile, assign=assignfile)

print("Inst. CG5")
cg = CG5(instance=inst, epslon=epslon,name=name)

print("build model")
# Initialize restricted master problem
cg.build_lp_model()

print("Loop")

# Initialize vars

k=0
impr=0

continue_cond = True

best_zrm = + np.inf
best_pi = np.zeros(inst.nproc)
best_alpha = np.zeros(inst.nmach)
best_omega = - np.inf

beta0 = .75
betaJ = .99
beta_min = 0
beta_step=10

beta = beta0

first_zrm = None

while continue_cond:
    start = time()
    # solve RDWM 
    print(" iter %5d         obj:   "% k ,end='',flush=True)
    (z_rm, pi_rm, alpha_rm, rlx ) = cg.solve_relax()

    if np.isnan(z_rm):
        print("relax unsolv")
        break

    if first_zrm is None:
        first_zrm = z_rm

    # Compute pi_stabilized
    pi = beta * best_pi + ( 1 - beta ) * pi_rm
    alpha = beta * best_alpha  + ( 1 - beta ) * alpha_rm

    # omega é a soma dos custos reduzidos dos sub problemas
    omega = 0
    print("%+17.3f (delta:  %17.3f O) in %10.3fs N: %3d, beta: %0.3f)" %( z_rm,  first_zrm - z_rm ,  rlx.Runtime, impr, beta) )

    if best_zrm > z_rm:
        best_zrm = z_rm 

    # solve L(pi_st)
    __col_added = False
    for m in range(inst.nmach):
        print(" iter %5d  mach %5d => " % (k,m),end='',flush=True)
        (roadef, w , q,model)= cg.compute_column(m, pi, alpha[m],k)
#        if w > - epslon:
#            print("for break")
##            break
        cg.validate_column(q,m,w,roadef,pi,alpha[m],epslon,model)
        omega += w
#        if w < -epslon or True:
#            __col_added = True
        s=cg.lp_add_col(obj=roadef, col =q, machine=m)
        print("%17.3f (roadef: %17.3f %s) in %10.3fs N: %3d, beta: %0.3f)" % (w, roadef,s, model.Runtime, impr, beta))
        if 'A' == s:
            __col_added = True

    print(" omega: %17.3f (%17.3f) delta: %17.3f  in %10.3fs N: %3d, beta: %0.3f)" % (omega, best_omega, omega - best_omega, time() - start , impr, beta))

    if omega > best_omega:
        impr += 1
        #beta = max(beta_min, beta0*(betaJ**max(0,impr - inst.nproc, k - inst.nproc)))
        print("     best omega improvement %17.3f -> %17.3f (%17.3f) %0.10f" % (best_omega, omega, omega/best_omega , beta ))
        #if np.isinf(best_omega):
        #    best_omega = omega
        #beta = beta * (omega/best_omega)
        best_omega = omega
        best_pi = pi
        best_alpha = alpha
        #input("Press Enter to continue...")
    #beta = max(beta_min, beta0*(betaJ**max(0,impr - inst.nmach)))
    #if impr > inst.nproc: beta = 0

    beta = max( beta_min, beta0 * (1 - (impr*1.0)/inst.nproc))

    #beta = beta0 * ((z_rm - omega)/(first_zrm - omega))

    if  best_omega > -epslon:
        contine_cond = False
        break
    if  not __col_added:
        print("  no cols added")
        contine_cond = False
        break

    
        
    k+=1
    if k>= 4 and False:
        break

(obj, X, alloc) = cg.solve_lp2mip()

print(X)
