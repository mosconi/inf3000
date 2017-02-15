#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from gurobipy import *

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-2)

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=True
epslon = 10^-1


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
    
    # Calculando custo da alocação inicial
    q[m].append(np.array([assign[p]==m for p in range(nproc)],dtype=np.int32))
    _util = (R*q[m][0].reshape(nproc,1)).sum(axis=0) 
    _obj1= _util - C_bar[m]
    _obj1[_obj1<0]=0

    _avail = C[m] - _util
    _obj2=np.array([[bT[r1,r2]*_avail[r1] - _avail[r2] for r2 in range(nres)] for r1 in range(nres)],dtype=np.int32)
    _obj2[_obj2<0] = 0
    
    _obj = (Wlc*_obj1).sum() + (Wbal*_obj2).sum()

    lbd[m].append(master_mdl.addVar(obj=_obj,vtype=GRB.BINARY,name="lbd_%d[0]"%m))

    
master_mdl.update()

# todo processo deve ser alocado
p_alloc=master_mdl.addConstrs(
    (quicksum(q[m][_a][p]*lbd[m][_a] for m in range(nmach) for _a in range(len(lbd[m]))) == 1 for p in range(nproc)), name="p_alloc")

# cada máquina somente pode possuir uma alocação
m_assign = master_mdl.addConstrs((quicksum(lbd[m][_a] for _a in range(len(lbd[m])))==1 for m in range(nmach)),name="m_assign")

mach_mdl={}
last_var=None
continue_condition = True
k = 0 
while continue_condition:
    
    master_mdl.update()     # atualiza antes de relaxar
    relax_mdl = master_mdl.relax()
    relax_mdl.Params.OutputFlag=0 # não imprime saída
    relax_mdl.optimize()

    rc = [(v.VarName, v.RC) for v in relax_mdl.getVars()]
    print("  RC", rc)
    pi = [relax_mdl.getConstrByName("p_alloc[%d]" % p ).Pi for p in range(nproc)]

    #print(pi[0])

    for m in range(nmach):
        print("maquina %d, round %d" % (m,k))
        
        alpha = relax_mdl.getConstrByName("m_assign[%d]"%m).Pi

        mach_mdl[m]=Model("machine_%d" % m)
        mach_mdl[m].ModelSense=GRB.MINIMIZE

        x=mach_mdl[m].addVars(nproc,vtype=GRB.BINARY,name="x")
        for p in range(nproc):
            x[p].start=0
            
        u=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="u",vtype=GRB.INTEGER)
        d=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="d",vtype=GRB.INTEGER)
        a=mach_mdl[m].addVars(nres,lb=0,ub=C[m],name="a",vtype=GRB.INTEGER)

        b=mach_mdl[m].addVars(nres,nres,lb=0,name="b",vtype=GRB.INTEGER)
        
        obj1 = mach_mdl[m].addVar(lb=0,name="obj1",vtype=GRB.INTEGER,obj=1)
        obj2 = mach_mdl[m].addVar(lb=0,name="obj2",vtype=GRB.INTEGER,obj=1)
        obj3 = mach_mdl[m].addVar(lb=0,name="obj3",vtype=GRB.INTEGER,obj=1)
        obj5 = mach_mdl[m].addVar(lb=0,name="obj5",vtype=GRB.INTEGER,obj=1)

        mach_mdl[m].update()


        mach_mdl[m].addConstrs((u[r] == quicksum(R[p,r]*x[p] for p in range(nproc)) for r in range(nres)), name="util_%d" % m)

        mach_mdl[m].addConstrs((a[r] == C[m,r] - u[r] for r in range(nres)), name="avail_%d" % m)

        mach_mdl[m].addConstrs((b[r1,r2] >= bT[r1,r2]*a[r1] - a[r2] for r1 in range(nres) for r2 in range(nres)), name="balance_%d" % m)
    
        mach_mdl[m].addConstrs((d[r] >= u[r] - C_bar[m,r] for r in range(nres)),name="overload_%d" % m)
        mach_mdl[m].update()

        mach_mdl[m].addConstr(obj1 == quicksum(Wlc[r]*d[r] for r in range(nres)))
        mach_mdl[m].addConstr(obj2 == quicksum(Wbal[r1,r2]*b[r1,r2] for r1 in range(nres) for r2 in range(nres)) )
        mach_mdl[m].addConstr(obj3 == quicksum(WPMC*RHO[p]*x[p] for p in [p for p in range(nproc) if p not in mach_assign[m]]))
        mach_mdl[m].addConstr(obj5 == quicksum(WMMC*MU[assign[p],m]*x[p] for p in [p for p in range(nproc) if p not in mach_assign[m]]))
        
        mach_mdl[m].setObjective(obj1 + obj2 + obj3 + obj5
                                 - quicksum((pi[p])*x[p] for p in [p for p in range(nproc) if p not in mach_assign[m]]) 
                                 + quicksum((WPMC*RHO[p])*(1-x[p]) for p in [p for p in range(nproc) if p in mach_assign[m]])
        )
        mach_mdl[m].Params.OutputFlag=0
        mach_mdl[m].optimize()
        #print("maq %d: obj: %0.2f alpha: %0.2f --> %0.2f" % (m,mach_mdl[m].ObjVal,alpha ,mach_mdl[m].ObjVal - alpha))
        #for r in range(nres):
            #if verbose: print("resource obj %d: %d" % (r, d[r].X))

        mach_mdl[m]._test = mach_mdl[m].ObjVal - alpha

        if mach_mdl[m].ObjVal > alpha:
            continue

        novo_q = np.array([1*(x[p].X>.5) for p in range(nproc)])
        if next((True for i in q[m] if np.array_equal(i, novo_q)), False):
            print("coluna já existe")
            mach_mdl[m]._skipped = True
            continue
        mach_mdl[m]._skipped = False

        _util = (R*novo_q.reshape(nproc,1)).sum(axis=0) 
        _obj1= _util - C_bar[m]
        _obj1[_obj1<0]=0
        
        _avail = C[m] - _util
        _obj2=np.array([[bT[r1,r2]*_avail[r1] - _avail[r2] for r2 in range(nres)] for r1 in range(nres)],dtype=np.int32)
        _obj2[_obj2<0] = 0
    
        _obj3 = RHO*novo_q

        _obj5 = MU[assign,m]

        _v1 = sum([pi[p]*novo_q[p] for p in range(nproc) if p not in mach_assign[m]])
        _v2 = sum([pi[p]*novo_q[p] for p in range(nproc)])
        _obj = (Wlc*_obj1).sum() + (Wbal*_obj2).sum() + (WPMC*_obj3).sum() + (WMMC*_obj5).sum()

        q[m].append(novo_q)
        col = Column()
   
        col.addTerms(q[m][-1],
                     [p_alloc[p] for p in range(nproc)])
        
        col.addTerms([1], [m_assign[m]])

        print("  obj calculado %12d"%_obj)
        print("  obj1          %12d"%(Wlc*_obj1).sum())
        print("  obj2          %12d"%(Wbal*_obj2).sum())
        print("  obj3          %12d"%(WPMC*_obj3).sum())
        print("  obj5          %12d"%(WMMC*_obj5).sum())
        print("  obj obtido       %12.2f"% (obj1.X + obj2.X + obj3.X + obj5.X))
        print("  obj1             %12.2f"% obj1.X)
        print("  obj2             %12.2f"% obj2.X)
        print("  obj3             %12.2f"% obj3.X)
        print("  obj5             %12.2f"% obj5.X)
        print("  v1               %12.2f"% _v1)
        print("  v2               %12.2f"% _v2)
        print("  alpha            %12.2f"% alpha)
        print("  delta:           %12.2f" % (_obj - (obj1.X + obj2.X + obj3.X + obj5.X)))
        lbd[m].append(master_mdl.addVar(obj=_obj,vtype=GRB.INTEGER,name="lbd_%d[%d]"%(m,len(lbd[m])), column=col))
        last_var=lbd[m][-1]
        # algo como : col.rc  <-- custo reduzido da coluna
        # col.rc == Obj(m)

    print([mach_mdl[m]._test for m in range(nmach)])

    continue_condition = all([mach_mdl[m]._test < - epslon for m in range(nmach)]) 
    k+=1

    if all([mach_mdl[m]._skipped for m in range(nmach)]):
        continue_condition = False



master_mdl.optimize()

for m in range(nmach):
    print("máquina %d allocs: %d"% (m, len(q[m])))
    print(np.array([[q[m][_a][p] for p in range(nproc)]for _a in range(len(q[m]))],dtype=np.int32))
    print(np.array([lbd[m][_a].X for _a in range(len(lbd[m]))]))
