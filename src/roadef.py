#!/usr/bin/env python3

import os, sys, getopt
import numpy as np

from gurobipy import *

def usage():
    print("%s [-n name] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"

try:
    opts, args = getopt.getopt(sys.argv[1:],"n:m:a:h",["name=","model=","assign=","help"])
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

for s in range(nserv):
    l = lines.pop(0)

nproc = lines.pop(0)[0]

S=[[] for i in range(nserv)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    PMC.append(l[0])

nbal = lines.pop(0)[0]
for b in range(nbal):
    l = lines.pop(0)
    l = lines.pop(0)
    
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
print(">>> machines ",nmach)
print(">>> process ",nproc)
print(">>> services ",nserv)
print(">>> neighborhood ",len(N))
print(">>> locations ",len(L))

mdl=Model(name)

mdl.ModelSense = GRB.MINIMIZE

print("creating x vars")
x = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,lb=0,ub=1,name="x")

print("creating z- vars")
z_minus = mdl.addVars(nproc,nmach,name="z-",vtype=GRB.INTEGER,lb=0,ub=1)

print("creating z+ vars")
z_plus = mdl.addVars(nproc,nmach,name="z+",vtype=GRB.INTEGER,lb=0,ub=1)

print("starting x[p,m], z-[p,m], z+[p,m]")
for p in range(nproc):
    for m in range(nmach):
        x[p,m].start = x_bar[p,m]
        z_minus[p,m].start =0
        z_plus[p,m].start =0
        
xi = mdl.addVars(nmach,nres,name="xi",lb=0)
u = mdl.addVars(nmach,nres,name="u",lb=0)

resobj = mdl.addVars(nres,name="res_obj",lb=0)
pmc = mdl.addVar(name="pmc",lb=0)


print("constr. all proc assigned")
mdl.addConstrs((x.sum(p,'*') == 1 for p in range(nproc)), name="process assigned")

print("constr. knapsack + xi")
one = [1 for m in range(nmach)]
for r in range(nres):
    for m in range(nmach):
        v = [x[p,m] for p in range(nproc)]
        mdl.addConstr(LinExpr(R[:,r],v) - u[m,r]==0, name="utilization")
        mdl.addConstr(u[m,r] <= C[m,r],name="cap")
        mdl.addConstr(u[m,r] - xi[m,r] <= SC[m,r],name="safecap")
    v = [xi[m,r] for m in range(nmach)]
    mdl.addConstr(resobj[r] == quicksum(v), name="obj")

print("constr. z*[p,m]")    

mdl.addConstrs((z_plus[p,m] >= x[p,m] - x_bar[p,m] for  p in range(nproc) for m in range(nmach) ))
mdl.addConstrs((z_minus[p,m] >= x_bar[p,m] - x[p,m] for  p in range(nproc) for m in range(nmach)) )

print("contr. service")    
for s in S:
    mdl.addConstrs((x.sum(p,'*') <=1 for p in s))


print("constr. PMC & MMC")
mdl.addConstr(pmc == quicksum(PMC[p]*z_plus[p,m] for m in range(nmach) for p in range(nproc)))

print("setting objs")

mdl.update()
mdl.setParam(GRB.Param.PoolSolutions, 100)
#mdl.Params.timeLimit=3600
mdl.NumObj=nres+1


for r in range(nres):
    mdl.setParam(GRB.Param.ObjNumber,r)
    mdl.ObjNPriority = 10
    mdl.ObjNWeight = Wlc[r]
    mdl.ObjNABSTol = sum(C[:,r])
    mdl.ObjNName = "resource " + str(r)
    resobj[r].ObjN=Wlc[r]

mdl.setParam(GRB.Param.ObjNumber,nres)
mdl.ObjNPriority = 1
#mdl.ObjNABSTol = 
mdl.ObjNName = "PMC"
mdl.ObjNWeight = WPMC
pmc.ObjN=WPMC

print("save model")
mdl.write(name + ".mps")
mdl.write(name + ".lp")

print("optimize")
mdl.optimize()

print('Optimization was stopped with status ' + str(mdl.Status))

for r in range(nres):
    print("resource obj %d: %d" % (r, resobj[r].X))
print("PMC            : %d" %pmc.X)

solution=[]
for p in range(nproc):
    for m in range(nmach):
        if x[p,m].X >0.5 :
            solution.append(m)

print(assign)
print(solution)

print([ "proc %d: %d -> %d" %(n,k[0],k[1])  for n,k in enumerate(zip(assign,solution)) if k[0]!=k[1]], sep="\n")

