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

x_bar= np.zeros((nproc,nmach), dtype=np.int32)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
assign=line[:]
for p in range(nproc):
    x_bar[p,line[p]] = 1

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

mdl_start = time()
load_time = mdl_start - all_begin

mdl=Model(name)
if not verbose: mdl.Params.OutputFlag=0

mdl.ModelSense = GRB.MINIMIZE

if verbose: print("creating x vars",flush=True)
x = mdl.addVars(nproc,nmach,vtype=GRB.BINARY,lb=0,ub=1,name="x")

if verbose: print("starting x[p,m] as x_bar[p,m]",flush=True)
for p in range(nproc):
    for m in range(nmach):
        x[p,m].start = x_bar[p,m]

if verbose: print("creating z- vars",flush=True)
z_minus = mdl.addVars(nproc,nmach,name="z-",vtype=GRB.BINARY)

for p in range(nproc):
    for m in range(nmach):
        z_minus[p,m].start=0

if verbose: print("creating t",flush=True)
t = mdl.addVars(nmach,nmach,vtype=GRB.SEMIINT,lb=0,name="t")

if verbose: print("creating u",flush=True)
u = mdl.addVars(nmach,nres,name="u",vtype=GRB.SEMIINT,lb=0,ub=C)

if verbose: print("creating a",flush=True)
a = mdl.addVars(nmach,nres,name="a",vtype=GRB.SEMIINT,lb=0,ub=C)

if verbose: print("creating d",flush=True)
d = mdl.addVars(nmach,nres,name="d",vtype=GRB.SEMIINT,lb=0)

if verbose: print("creating b",flush=True)
b = mdl.addVars(nmach,nres,nres,vtype=GRB.SEMIINT,name="b")

if verbose: print("creating o",flush=True)
o = mdl.addVars(nserv,len(L),vtype=GRB.BINARY,name="o")

if verbose: print("creating g",flush=True)
g = mdl.addVars(nserv,vtype=GRB.BINARY,name="s")

if verbose: print("creating h",flush=True)
h = mdl.addVars(nserv,len(N),vtype=GRB.BINARY,name="h")

for s in range(nserv):
    for n in range(len(N)):
        h[s,n].start=0

if verbose: print("creating obj vars",flush=True)
loadcost = mdl.addVars(nres,name="locadcost",vtype=GRB.SEMIINT,obj=Wlc)
balancecost= mdl.addVars(nres,nres,name="balancecost",vtype=GRB.SEMIINT,obj=Wbal)
pmc = mdl.addVar(name="pmc",vtype=GRB.SEMIINT,obj=WPMC)
smc = mdl.addVar(name="smc",vtype=GRB.SEMIINT,obj=WSMC)
mmc = mdl.addVar(name="mmc",vtype=GRB.SEMIINT,obj=WMMC)

mdl_vars = time() - mdl_start

if verbose: print("model update",flush=True)
mdl.update()

if verbose: print("constr. all proc assigned")
mdl.addConstrs((x.sum(p,'*') == 1 for p in range(nproc)), name="process_assigned")

if verbose: print("constr. utilization")
mdl.addConstrs((u[m,r] == quicksum(x[p,m]*R[p,r] for p in range(nproc)) for r in range(nres) for m in range(nmach)),name="utilization")

if verbose: print("constr. util < cap")
mdl.addConstrs((u[m,r] <= C[m,r] for r in range(nres) for m in range(nmach)),
               name="cap")

if verbose: print("constr. avail")
mdl.addConstrs((a[m,r] == C[m,r] - u[m,r] for r in range(nres) for m in range(nmach)),
               name="avail")

if verbose: print("constr. transient")
for r in range(nres):
    if T[r]:
        mdl.addConstrs((quicksum(x[p,m]*R[p,r] for p in range(nproc)) + quicksum(z_minus[p,m]*R[p,r] for p in range(nproc)) <= C[m,r] for m in range(nmach)),name="transient")

if verbose: print("constr. overload")
mdl.addConstrs((u[m,r] - d[m,r] <= C_bar[m,r] for r in range(nres) for m in range(nmach)),name="overload")


if verbose: print("constr. z[p,m]")    

mdl.addConstrs((z_minus[p,m] >= x_bar[p,m] - x[p,m] for  p in range(nproc) for m in range(nmach)) )

if verbose: print("constr. conflito")
mdl.addConstrs((quicksum(x[p,m] for p in S[s])<=1 for s in range(len(S)) for m in range(nmach)), name="conflito")

if verbose: print("constr. t[i,j]=sum(x_bar[p,i]*x[p,j])")
mdl.addConstrs((t[i,j]==quicksum(x_bar[p,i]*x[p,j] for p in range(nproc)) for i in range(nmach) for j in range(nmach)),name="t")

if verbose: print("constr. o[s,l]")
for s in range(nserv):
    for l in range(len(L)):
        mdl.addConstrs((o[s,l] >= x[p,m] for p in S[s] for m in L[l]),name=("o[%d,%d]"%(s,l)))
        mdl.addConstr((o[s,l]<=quicksum((x[p,m] for p in S[s] for m in L[l]))), name="o_ub")
        
mdl.addConstrs((quicksum(o[s,l] for l in range(len(L))) >= delta[s] for s in range(nserv)),name="serv_spread")

if verbose: print("constr. g[s]")
mdl.addConstrs((g[s] == quicksum(z_minus[p,m] for m in range(nmach) for p in S[s]) for s in range(nserv)), name="g")

if verbose: print("constr. h[s,n]")
for s in range(nserv):
    for n in range(len(N)):
        mdl.addConstrs((h[s,n] >= x[p,m] for p in S[s] for m in N[n]),name=("h[%d,%d]"%(s,l)))
        mdl.addConstr((h[s,n] <= quicksum((x[p,m] for p in S[s] for m in N[n]))),name="h_upperbound")

for s in range(nserv):
    if s in sdep:
        for _s in sdep[s]:
            mdl.addConstrs((h[s,n] <= h[_s,n] for n in range(len(N))),name=("dep[%d,%d]"%(s,_s)))

if verbose: print("constr. b[m,r1,r2]")
mdl.addConstrs((b[m,r1,r2] >= bT[r1,r2]*a[m,r1] - a[m,r2] for m in range(nmach) for r1 in range(nres) for r2 in range(nres)), name="b")    


if verbose: print("constr. loadcost")
mdl.addConstrs((loadcost[r] == quicksum(d[m,r] for m in range(nmach)) for r in range(nres)), name="loadcost")

mdl.addConstrs((balancecost[r1,r2] == quicksum(b[m,r1,r2] for m in range(nmach)) for r1 in range(nres) for r2 in range(nres)))

if verbose: print("constr. PMC")
mdl.addConstr(pmc == quicksum(RHO[p]*z_minus[p,m] for m in range(nmach) for p in range(nproc)))

if verbose: print("constr. SMC")
mdl.addConstrs(smc >= g[s] for s in range(nserv))

if verbose: print("constr. MMC")
mdl.addConstr(mmc == quicksum(MU[i,j]*t[i,j] for i in range(nmach) for j in range(nmach)))

mdl_constrs = time() - mdl_start - mdl_vars

if savemodel:
    if verbose: print("save model")
    mdl.write(name + ".mps")
    mdl.write(name + ".lp")

def cb(model,where):
    """ Callback

    O modelo terminará antes, aceitando um gap de 5%, pois queremos medir a 
    qualidade do modelo, em relação a fornecer soluções viáveis.  
    Posteriormente será removido este callback.
    """
    if where == GRB.Callback.MIP:
        # General MIP callback
        objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
        objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
        if abs(objbst - objbnd) < 0.05 * (1.0 + abs(objbst)):
            if verbose: print('>>>> Stop early - 5% gap achieved')
            model.terminate()

mdl._x=x
if verbose: print("optimize")
_start = time()
mdl.optimize(cb)
_elapsed = time() - _start
mdl_opt=_elapsed

if verbose: print('Optimization was stopped with status ' + str(mdl.Status),flush=True)

if verbose: print(">>> optimized in %0.2f" % _elapsed ,flush=True)

if verbose: print(">>> Optimal value: %d" % mdl.objVal)

for r in range(nres):
    if verbose: print("resource obj %d: %d" % (r, loadcost[r].X))
for r1 in range(nres):
    for r2 in range(nres):
        if verbose: print("balance (%d,%d): %d"%(r1,r2,balancecost[r1,r2].X))
if verbose: print("PMC            : %d" %pmc.X)
if verbose: print("SMC            : %d" %smc.X)
if verbose: print("MMC            : %d" %mmc.X)

solution=[]
for p in range(nproc):
    for m in range(nmach):
        if x[p,m].X >0.5 :
            solution.append(m)

if verbose: print(assign)
if verbose: print(solution)

if verbose: print([ "proc %d: %d -> %d" %(n,k[0],k[1])  for n,k in enumerate(zip(assign,solution)) if k[0]!=k[1]], sep="\n")

if outputfile:
    if verbose: print(">>> salvando resposta em %s" % outputfile)
    f = open(outputfile, "w")
    f.write(' '.join(str(i) for i in solution))
    f.close()

all_time = time() - all_begin
if tex:
    import subprocess
    import resource
    output = subprocess.check_output(["./checker", modelfile,assignfile,outputfile])
    print("\\verb|%s| & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %d & %0.3f & %d & %0.3f  & %s \\\\" %
          (name, all_time, load_time, mdl_vars, mdl_constrs, mdl_opt,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 +1 , mdl._objbst, mdl.objVal,  (abs(mdl._objbst-mdl.objVal)/abs(mdl._objbst))*100 ,output.decode("utf-8").strip('\n')))
