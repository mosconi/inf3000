#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from gurobipy import *

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=False

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

all_begin = time()

T=[]
Wlc=[]
M=[]
C=[]
C_bar=[]
S=[]
R=[]
RHO=[]
MU=[]

f = open(modelfile, "r")

lines = [list(map(int,l.rstrip('\n').split())) for l in f]

nres = lines.pop(0)[0]

for r in range(nres):
    l = lines.pop(0)
    T.append(l[0])
    Wlc.append(l[1])


Wlc=np.array(Wlc)
nmach = lines.pop(0)[0]

L=[[] for i in range(nmach)]
N=[[] for i in range(nmach)]

for m in range(nmach):
    l = lines.pop(0)
    N[l.pop(0)].append(m) # neighborhood
    L[l.pop(0)].append(m) # location
    C.append(l[:nres])
    del(l[:nres])
    C_bar.append(l[:nres])
    del(l[:nres])
    MU.append(l)

# remove empty neighborhood/locations
L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
C_bar=np.array(C_bar)
MU=np.array(MU)

nserv = lines.pop(0)[0]
delta=[0 for s in range(nserv)]
sdep={}

for s in range(nserv):
    l = lines.pop(0)
    delta[s]=l[0]
    if l[1]>0:
        sdep[s]=l[2:]

nproc = lines.pop(0)[0]

S=[[] for i in range(nserv)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    RHO.append(l[0])

nbal = lines.pop(0)[0]
bT=np.zeros((nres,nres),dtype=np.int32)
Wbal=np.zeros((nres,nres),dtype=np.int32)
for b in range(nbal):
    l = lines.pop(0)
    r1=l[0]
    r2=l[1]
    bT[ r1,r2]=l[2]
    l = lines.pop(0)
    Wbal[r1,r2] = l[0]

WPMC=lines[0][0]
WSMC=lines[0][1]
WMMC=lines[0][2]
    
R = np.array(R)
S = [s for s in S if s]    

x0= np.zeros((nproc,nmach), dtype=np.int32)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
assign=line[:]
for p in range(nproc):
    x0[p,line[p]] = 1

T=np.array(T)

if verbose: print(">> problem size")
if verbose: print(">>> resources ",nres)
if verbose: print(">>> transisent ",T.sum())
if verbose: print(">>> machines ",nmach)
if verbose: print(">>> process ",nproc)
if verbose: print(">>> services ",nserv)
if verbose: print(">>> balances ",nbal)
if verbose: print(">>> neighborhood ",len(N))
if verbose: print(">>> locations ",len(L))

load_time = time() - all_begin

mdl=Model(name)
if not verbose: mdl.Params.OutputFlag=0

mdl.ModelSense = GRB.MINIMIZE

mdl_start = time()
_start = time()

if verbose: print("creating x vars",flush=True)
x = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,lb=0,ub=1,name="x")

if verbose: print("starting x[p,m] as x0[p,m]",flush=True)
for p in range(nproc):
    for m in range(nmach):
        x[p,m].start = x0[p,m]

if verbose: print("creating z vars",flush=True)
z = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,lb=0,ub=1,name="z")

if verbose: print("creating u",flush=True)
u = mdl.addVars(nmach,nres,name="u",vtype=GRB.INTEGER,lb=0,ub=C)
ut = mdl.addVars(nmach,nres,name="ut",vtype=GRB.INTEGER,lb=0,ub=C)

if verbose: print("creating d",flush=True)
d = mdl.addVars(nmach,nres,name="d",vtype=GRB.INTEGER,lb=0,ub=C)

mdl.update()

if verbose: print("constr. utilization")
mdl.addConstrs((u[m,r] == quicksum(x[p,m]*R[p,r] for p in range(nproc)) for r in range(nres) for m in range(nmach)),name="utilization")

if verbose: print("constr. z")
mdl.addConstrs((z[p,m] >= x0[p,m] - x[p,m] for p in range(nproc) for m in range(nmach)),name="z")

if verbose: print("constr. utilization transient")
mdl.addConstrs((ut[m,r] == 
                quicksum(z[p,m]*R[p,r]*T[r] for p in range(nproc))
                for r in range(nres) for m in range(nmach)),name="utilization_transient")

if verbose: print("constr. all proc assigned")
mdl.addConstrs((quicksum(x[p,m] for m in range(nmach)) == 1 for p in range(nproc)), name="process_assigned")


if verbose: print("constr. util + ut <= cap (hard) + (hard)")
mdl.addConstrs((u[m,r] + ut[m,r] <= C[m,r] for r in range(nres) for m in range(nmach)),
               name="cap")

if verbose: print("constr. overload (obj)")
mdl.addConstrs((u[m,r] - d[m,r] <= C_bar[m,r] for r in range(nres) for m in range(nmach)),name="overload")


mdl.setObjective(quicksum(Wlc[r]*d[m,r] for m in range(nmach) for r in range(nres)) +
                 quicksum(RHO[p]*z[p,m] for m in range(nmach) for p in range(nproc)) , GRB.MINIMIZE)

mdl.update()
mdl.optimize()

def validate_solution(X,model):
    for s in (s for s in range(nserv) if len(S[s])>1):
        for m in range(nmach):
            if sum(x[p,m].X for p in S[s])>1:
                print("conflito do serviço %s falhou em %d" % (s,m))
        o={}
        for l in range(len(L)):
            o[s,l] = any(x[p,m] for p in S[s] for m in L[l])
            print("o[%d,%d] = %d" % (s,l,o[s,l]))
        print("sum(o[%d,*]) = %d" % (s,sum(o[s,l] for l in range(len(L)))))
        if sum(o[s,l] for l in range(len(L))) < delta[s]:
            print("spread do serviço %s falhou %d < %d" %(s,sum(o[s,l] for l in range(len(L))), delta[s]))
            

validate_solution(x,mdl)

print(" obj = %0.3f" % mdl.objVal)

solution=[]
for p in range(nproc):
    for m in range(nmach):
        if x[p,m].X >0.5 :
            solution.append(m)

if verbose: print(assign)
if verbose: print(solution)

if outputfile:
    if verbose: print(">>> salvando resposta em %s" % outputfile)
    f = open(outputfile, "w")
    f.write(' '.join(str(i) for i in solution))
    f.close()

if verbose: print([ "proc %d: %d -> %d" %(n,k[0],k[1])  for n,k in enumerate(zip(assign,solution)) if k[0]!=k[1]], sep="\n")
