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


best_omega =  -float('Inf')

continue_cond = True

best_pi = np.zeros(inst.nproc)
best_alpha = np.zeros(inst.nmach)
best_omega = - np.inf

beta = .5

cg.build_model()

k=0

while continue_cond:
    (z_rm, pi_rm, alpha_rm) = cg.solve_relax()

    pi = beta * best_pi + ( 1 - beta ) * pi_rm
    alpha = beta * best_alpha  + ( 1 - beta ) * alpha_rm

    omega = 0
    print(" iter %05d  obj: %+20.3f" %(k, z_rm) )
    
    for m in range(inst.nmach):
        print(" iter %05d  mach %05d => " % (k,m),end='')
        (roadef, w , q)= cg.compute_column(m, pi, alpha[m])
#        if w > - epslon:
#            print("for break")
#            break

        print("%20.3f (roadef: %20.3f rat: %20.10f)" % (w, roadef,abs(w/(roadef+1))))
        omega += w
        if len(q) < 100:
            print("add col", q)
        cg.mip_add_col(obj=roadef, col =q, machine=m)

        #cols = cg.generate_companion_columns(machine=m,col=q,obj=roadef
        #for (_m,_obj,_col) in cols 

    print("omega: %20.3f  (best:  %20.3f delta: %20.3f rat(obj/best): %20.3f)" % (omega, best_omega, z_rm - best_omega, z_rm/best_omega))

    if omega > best_omega:
        print("best omega improvement %20.3f -> %20.3f" % (best_omega, omega))
        best_omega = omega
        best_pi = pi_rm
        best_alpha = alpha_rm
        input("Press Enter to continue...")

    if abs(z_rm - best_omega) < epslon:
        contine_cond = False
        break
        
    k+=1


(obj, X, alloc) = cg.solve_mip()

print(X)
