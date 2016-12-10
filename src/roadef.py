#!/usr/bin/env python3

import os, sys, getopt
import numpy as np

from gurobipy import *

def usage():
    print("%s -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None


try:
    opts, args = getopt.getopt(sys.argv[1:],"m:a:h",["model=","assign=","help"])
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

L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
SC=np.array(SC)
MMC=np.array(M)

nserv = lines.pop(0)[0]

for s in range(nserv):
    l = lines.pop(0)

nproc = lines.pop(0)[0]

S=[[] for i in range(nproc)]

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

mdl=Model("roadef")

mdl.ModelSense = GRB.MINIMIZE

print("creating x vars")
x = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,name="x")

print("creating z- vars")
z_minus = mdl.addVars(nproc,nmach,name="z-",vtype=GRB.INTEGER)

print("creating z+ vars")
z_plus = mdl.addVars(nproc,nmach,name="z+",vtype=GRB.INTEGER)

#print("creating mig vars")
#mig = mdl.addVars(nproc,nmach,nmach,name="mig",vtype=GRB.INTEGER)

#print("creating z vars")
#z = mdl.addVars(nproc,name="z",vtype=GRB.BINARY)

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
#mmc = mdl.addVar(name="mmc",lb=0)

mdl.update()

print("constr. all proc assigned")
mdl.addConstrs((x.sum(p,'*') == 1 for p in range(nproc)), name="all proc assigned")

print("knapsack")
one = [1 for m in range(nmach)]
for r in range(nres):
    for m in range(nmach):
        v = [x[p,m] for p in range(nproc)]
        mdl.addConstr(LinExpr(R[:,r],v) == u[m,r])
        mdl.addConstr(u[m,r] <= C[m,r])
        mdl.addConstr(u[m,r] - xi[m,r] <= SC[m,r])
    v = [xi[m,r] for m in range(nmach)]
    mdl.addConstr(resobj[r] == LinExpr(one,v))

print("z*[p,m]")    
for p in range(nproc):
    for m in range(nmach):
        mdl.addConstr(z_plus[p,m] >= x[p,m] - x_bar[p,m])
        mdl.addConstr(z_minus[p,m] >= x_bar[p,m] - x[p,m])

print("service constraint")    
for s in S:
    mdl.addConstrs((x.sum(p,'*') <=1 for p in s))

#print("migration")
#for i in range(nmach):
#    for j in range(nmach):
#        for p in range(nproc):
#            mdl.addConstr(mig[p,i,j] >= x_bar[p,i] +  x[p,j] - 1)
#            mdl.addConstr(mig[p,i,j] <= x_bar[p,i])
#            mdl.addConstr(mig[p,i,j] <= x[p,j])

print("PMC & MMC")
mdl.addConstr(pmc == quicksum(PMC[p]*z_plus[p,m] for m in range(nmach) for p in range(nproc)))
#mdl.addConstr(mmc == quicksum(MMC[i,j]*mig[p,i,j] for i in range(nmach) for j in range(nmach) for p in range(nproc)))


print("setting objs")

mdl.update()
mdl.setParam(GRB.Param.PoolSolutions, 100)
#mdl.Params.timeLimit=300
mdl.NumObj=nres+1


for r in range(nres):
    mdl.setParam(GRB.Param.ObjNumber,r)
    mdl.ObjNPriority = 10
#    mdl.ObjNWeight = Wlc[r]
    mdl.ObjNRelTol = 1
    mdl.ObjNName = "resource " + str(r)
    resobj[r].ObjN=1

mdl.setParam(GRB.Param.ObjNumber,nres)
mdl.ObjNPriority = 1
mdl.ObjNRelTol = 1
mdl.ObjNName = "PMC"
#mdl.ObjNWeight = WPMC
pmc.ObjN=1

#mdl.setParam(GRB.Param.ObjNumber,nres+1)
#mdl.ObjNPriority = 1
#mdl.ObjNRelTol = 1
#mdl.ObjNName = "MMC"
#mdl.ObjNWeight = WMMC
#mmc.ObjN=1

print("optimize")
mdl.optimize()

print('Optimization was stopped with status ' + str(mdl.Status))

for r in range(nres):
    print("resource obj %d: %d" % (r, resobj[r].X))
print("PMC            : %d" %pmc.X)
#print("MMC            : %d" %mmc.X)

solution=[]
for p in range(nproc):
    for m in range(nmach):
        if x[p,m].X ==1:
            solution.append(m)

print(assign)
print(solution)
