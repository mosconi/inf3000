#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from gurobipy import *

def usage():
    print("%s [-n name] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None

try:
    opts, args = getopt.getopt(sys.argv[1:],"n:m:a:o:h",["name=","model=","assign=","output=","help"])
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

T=[]
Wlc=[]
M=[]
C=[]
SC=[]
S=[]
R=[]
PMC=[]
MMC=[]

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
    SC.append(l[:nres])
    del(l[:nres])
    M.append(l)

# remove empty neighborhood/locations
L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
SC=np.array(SC)
MMC=np.array(M)

nserv = lines.pop(0)[0]
spread=[0 for s in range(nserv)]
sdep={}

for s in range(nserv):
    l = lines.pop(0)
    spread[s]=l[0]
    if l[1]>0:
        sdep[s]=l[2:]

nproc = lines.pop(0)[0]

S=[[] for i in range(nserv)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    PMC.append(l[0])

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

x_bar= np.zeros((nproc,nmach), dtype=np.int32)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
assign=line[:]
for p in range(nproc):
    x_bar[p,line[p]] = 1

T=np.array(T)

print(">> problem size")
print(">>> resources ",nres)
print(">>> transisent ",T.sum())
print(">>> machines ",nmach)
print(">>> process ",nproc)
print(">>> services ",nserv)
print(">>> balances ",nbal)
print(">>> neighborhood ",len(N))
print(">>> locations ",len(L))

mdl=Model(name)

mdl.ModelSense = GRB.MINIMIZE

print("creating x vars",flush=True)
x = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,lb=0,ub=1,name="x")

print("starting x[p,m] as x_bar[p,m]",flush=True)
for p in range(nproc):
    for m in range(nmach):
        x[p,m].start = x_bar[p,m]

print("creating z- vars",flush=True)
z_minus = mdl.addVars(nproc,nmach,name="z-",vtype=GRB.INTEGER,lb=0,ub=1)

print("creating z+ vars",flush=True)
z_plus = mdl.addVars(nproc,nmach,name="z+",vtype=GRB.INTEGER,lb=0,ub=1)


print("creating y",flush=True)
y = mdl.addVars(nproc,nmach,nmach,vtype=GRB.INTEGER,lb=0,ub=1,name="y")

print("creating t",flush=True)
t = mdl.addVars(nmach,nmach,vtype=GRB.SEMIINT,lb=0,name="t")

print("creating o",flush=True)
o = mdl.addVars(nserv,len(L),vtype=GRB.BINARY,name="o")

print("creating g",flush=True)
g = mdl.addVars(nserv,vtype=GRB.BINARY,name="s")

print("creating h",flush=True)
h = mdl.addVars(nserv,len(N),vtype=GRB.BINARY,name="h")

print("creating k",flush=True)
k = mdl.addVars(nproc,len(N),vtype=GRB.BINARY,name="k")

print("creating b",flush=True)
b = mdl.addVars(nmach,nres,nres,vtype=GRB.SEMIINT,name="b")

print("creating d",flush=True)
d = mdl.addVars(nmach,nres,name="d",vtype=GRB.SEMIINT,lb=0)

print("creating u",flush=True)
u = mdl.addVars(nmach,nres,name="u",vtype=GRB.SEMIINT,lb=0,ub=C)

print("creating a",flush=True)
a = mdl.addVars(nmach,nres,name="a",vtype=GRB.SEMIINT,lb=0,ub=C)


print("creating obj vars",flush=True)
loadcost = mdl.addVars(nres,name="locadcost",vtype=GRB.SEMIINT,obj=Wlc)
balancecost= mdl.addVars(nres,nres,name="balancecost",vtype=GRB.SEMIINT,obj=Wbal)
pmc = mdl.addVar(name="pmc",vtype=GRB.SEMIINT,obj=WPMC)
smc = mdl.addVar(name="smc",vtype=GRB.SEMIINT,obj=WSMC)
mmc = mdl.addVar(name="mmc",vtype=GRB.SEMIINT,obj=WMMC)

print("model update",flush=True)
mdl.update()


print("constr. all proc assigned")
mdl.addConstrs((x.sum(p,'*') == 1 for p in range(nproc)), name="process_assigned")

print("constr. utilization")
mdl.addConstrs((u[m,r] == quicksum(x[p,m]*R[p,r] for p in range(nproc)) for r in range(nres) for m in range(nmach)),name="utilization")

print("constr. util < cap")
mdl.addConstrs((u[m,r] <= C[m,r] for r in range(nres) for m in range(nmach)),
               name="cap")

print("constr. avail")
mdl.addConstrs((a[m,r] == C[m,r] - u[m,r] for r in range(nres) for m in range(nmach)),
               name="avail")

print("constr. transient")
for r in range(nres):
    if T[r]:
        mdl.addConstrs((quicksum(x[p,m]*R[p,r] for p in range(nproc)) + quicksum(z_minus[p,m]*R[p,r] for p in range(nproc)) <= C[m,r] for m in range(nmach)),name="transient")

print("constr. overload")
mdl.addConstrs((u[m,r] - d[m,r] <= SC[m,r] for r in range(nres) for m in range(nmach)),name="overload")

print("constr. loadcost")
mdl.addConstrs((loadcost[r] == quicksum(d[m,r] for m in range(nmach)) for r in range(nres)), name="loadcost")

print("constr. z*[p,m]")    

mdl.addConstrs((z_plus[p,m] >= x[p,m] - x_bar[p,m] for  p in range(nproc) for m in range(nmach) ))
mdl.addConstrs((z_minus[p,m] >= x_bar[p,m] - x[p,m] for  p in range(nproc) for m in range(nmach)) )

print("constr. services")
mdl.addConstrs((quicksum(x[p,m] for p in S[s])<=1 for s in range(len(S)) for m in range(nmach)), name="service")

print("constr. y[p,i,j] >= z- + z+")
mdl.addConstrs((y[p,i,j] >= z_minus[p,i] + z_plus[p,j] - 1 for p in range(nproc) for i in range(nmach) for j in range(nmach)), name="y")
#print("constr. y[p,i,j] <= z-")
#mdl.addConstrs((y[p,i,j] <= z_minus[p,i] for p in range(nproc) for i in range(nmach) for j in range(nmach)), name="y_z-")
#print("constr. y[p,i,j] <= z+")
#mdl.addConstrs((y[p,i,j] <= z_plus[p,j]  for p in range(nproc) for i in range(nmach) for j in range(nmach)), name="y_z+")

print("constr. t[i,j]=sum(y[p,i,j])")
mdl.addConstrs((t[i,j]==quicksum(y[p,i,j] for p in range(nproc)) for i in range(nmach) for j in range(nmach)),name="t")

print("constr. o[s,l]")
for s in range(nserv):
    for l in range(len(L)):
        mdl.addConstrs((o[s,l] >= x[p,m] for p in S[s] for m in L[l]),name=("o[%d,%d]"%(s,l)))

mdl.addConstrs((quicksum(o[s,l] for l in range(len(L))) >= spread[s]) for s in range(nserv))

#print("constr. k[p,n]")
#mdl.addConstrs((k[p,n] == quicksum(x[p,m] for m in N[n]) for p in range(nproc) for n in range(len(N))),name="k")

print("constr. g[s]")
mdl.addConstrs((g[s] == quicksum(z_plus[p,m] for m in range(nmach) for p in S[s]) for s in range(nserv)), name="g")

print("constr. h[s,n]")
for s in range(nserv):
    for n in range(len(N)):
        mdl.addConstrs((h[s,n] >= x[p,m] for p in S[s] for m in N[n]),name=("h[%d,%d]"%(s,l)))

for s in range(nserv):
    if s in sdep:
        print((s,sdep[s]))
        for _s in sdep[s]:
            mdl.addConstrs((h[s,n] <= h[_s,n] for n in range(len(N))),name=("dep[%d,%d]"%(s,_s)))

print("constr. b[m,r1,r2]")
mdl.addConstrs((b[m,r1,r2] >= bT[r1,r2]*a[m,r1] - a[m,r2] for m in range(nmach) for r1 in range(nres) for r2 in range(nres)), name="b")    

mdl.addConstrs((balancecost[r1,r2] == quicksum(b[m,r1,r2] for m in range(nmach)) for r1 in range(nres) for r2 in range(nres)))

print("constr. PMC")
mdl.addConstr(pmc == quicksum(PMC[p]*z_plus[p,m] for m in range(nmach) for p in range(nproc)))

print("constr. SMC")
mdl.addConstrs(smc >= g[s] for s in range(nserv))

print("constr. MMC")
mdl.addConstr(mmc == quicksum(MMC[i,j]*t[i,j] for i in range(nmach) for j in range(nmach)))


#mdl.Params.OutputFlag=0
#mdl.Params.PoolSolutions = 1000
#mdl.Params.Threads = 1
#mdl.Params.MIPFocus = 3
#mdl.Params.Heuristics = 1

print("save model")
mdl.write(name + ".mps")
mdl.write(name + ".lp")
mdl._lastiter = -GRB.INFINITY
mdl._lastnode = -GRB.INFINITY


def cb(model,where):
    if where == GRB.Callback.MIP:
        # General MIP callback
        nodecnt = model.cbGet(GRB.Callback.MIP_NODCNT)
        objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
        objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
        solcnt = model.cbGet(GRB.Callback.MIP_SOLCNT)
        if abs(objbst - objbnd) < 0.05 * (1.0 + abs(objbst)):
            print('>>>> Stop early - 5% gap achieved')
            model.terminate()

mdl._x=x
print("optimize")
_start = time()
mdl.optimize(cb)
_elapsed = time() - _start

print('Optimization was stopped with status ' + str(mdl.Status),flush=True)

print(">>> optimized in %0.2f" % _elapsed ,flush=True)

for r in range(nres):
    print("resource obj %d: %d" % (r, loadcost[r].X))
for r1 in range(nres):
    for r2 in range(nres):
        print("balance (%d,%d): %d"%(r1,r2,balancecost[r1,r2].X))
print("PMC            : %d" %pmc.X)
print("SMC            : %d" %smc.X)
print("MMC            : %d" %mmc.X)

solution=[]
for p in range(nproc):
    for m in range(nmach):
        if x[p,m].X >0.5 :
            solution.append(m)

print(assign)
print(solution)

print([ "proc %d: %d -> %d" %(n,k[0],k[1])  for n,k in enumerate(zip(assign,solution)) if k[0]!=k[1]], sep="\n")

if outputfile:
    print(">>> salvando resposta em %s" % outputfile)
    f = open(outputfile, "w")
    f.write(' '.join(str(i) for i in solution))
