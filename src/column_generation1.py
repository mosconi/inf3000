#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from gurobipy import *

np.set_printoptions(linewidth=200)

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=True

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

RHO=np.array(RHO)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
assign=line[:]

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


# converte o vetor de alocação idexado por p, em listas de processos em M

mach_assign = [[p[0] for p in enumerate(assign) if p[1] == m ] for m in range(nmach)]

master_mdl=Model("master")
master_mdl.ModelSense=GRB.MINIMIZE

lbd=[[] for m in range(nmach)]
q=[[] for m in range(nmach)]

for m in range(nmach):
    lbd[m].append(np.array([assign[p]==m for p in range(nproc)],dtype=np.int32))
    q[m].append(master_mdl.addVars(nproc,obj=1,vtype=GRB.BINARY,name="q_%d[0]"%m))

master_mdl.update()
c_alloc=master_mdl.addConstrs(
    (quicksum(q[m][_a][p]*lbd[m][_a][p] for m in range(nmach) for _a in range(len(q[m]))) == 1 for p in range(nproc)),
    name="alloc")

master_mdl.update()

mach_mdl={}
for k in range(5):

    relax_mdl = master_mdl.relax()

    relax_mdl.Params.OutputFlag=0
    relax_mdl.optimize()

    pi = [c.Pi for c in relax_mdl.getConstrs()]

    print(len(pi))
    print(pi)

    for m in range(nmach):
        print("maquina %d" % m)
        mach_mdl[m]=Model("machine_%d" % m)
        mach_mdl[m].ModelSense=GRB.MINIMIZE

        x=mach_mdl[m].addVars(nproc,vtype=GRB.BINARY,name="x")
        for p in range(nproc):
            x[p].start=0
            
        u=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="u",vtype=GRB.INTEGER)
        d=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="d",vtype=GRB.INTEGER,obj=Wlc)
        a=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="a",vtype=GRB.INTEGER)

        b=mach_mdl[m].addVars(nres,nres,lb=0,name="b",vtype=GRB.INTEGER,obj=Wbal)

        mach_mdl[m].update()

        print([ sum(R[p,r]*x[p].start for p in range(nproc)) -C_bar[m,r] for r in range(nres)])

        beta= np.array(RHO,copy=True) 
        for p in [i for i in range(nproc) if i not in mach_assign[m]]:
            beta[p]+=MU[assign[p],m] - pi[p]

        mach_mdl[m].addConstrs((u[r] == quicksum(R[p,r]*x[p] for p in range(nproc)) for r in range(nres)), name="util_%d" % m)

        mach_mdl[m].addConstrs((a[r] == C[m,r] - u[r] for r in range(nres)), name="avail_%d" % m)

        mach_mdl[m].addConstrs((b[r1,r2] >= bT[r1,r2]*a[r1] - a[r2] for r1 in range(nres) for r2 in range(nres)), name="balance_%d" % m)
    
        mach_mdl[m].addConstrs((d[r] >= u[r] - C_bar[m,r] for r in range(nres)),name="overload_%d" % m)
        mach_mdl[m].update()
    
        mach_mdl[m].setObjective(
            quicksum(d[r] for r in range(nres)) +
            quicksum(b[r1,r2] for r1 in range(nres) for r2 in range(nres)) +
            quicksum(beta[p]*x[p] for p in [p for p in range(nproc) if p not in mach_assign[m]]) +
            1
        )
        mach_mdl[m].Params.OutputFlag=0
        mach_mdl[m].optimize()
        print("maq %d: obj: %0.2f" % (m,mach_mdl[m].ObjVal))
        for r in range(nres):
            if verbose: print("resource obj %d: %d" % (r, d[r].X))

        lbd[m].append(np.array([1*(x[p].X>.5) for p in range(nproc)]))
        col = Column()
        col.addTerms(lbd[m][-1],
                     [c_alloc[p] for p in range(nproc)])
        q[m].append(master_mdl.addVar(obj=1,vtype=GRB.INTEGER,name="q_%d[%d"%(m,len(lbd)-1)))
        master_mdl.update()

    master_mdl.optimize()

for m in range(nmach):
    print("máquina %d"%m)
    print(len(lbd[m]))
    for _a in range(len(lbd[m])):
        print(lbd[m][_a])
