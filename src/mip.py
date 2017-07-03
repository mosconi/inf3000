import argparse

import roadef
import common
import numpy as np

from gurobipy import *

class MipArgs(object):
    pass

args = MipArgs()


parser = argparse.ArgumentParser(add_help=False,
    description="Arguments required for model/solver")

solver = parser.add_argument_group("solver","Options to change solver behaviour")
solver.add_argument("--epslon",type=float,default=10**-1,
                    help="Error tolerance")
solver.add_argument("--runname",dest="run_name",default=None,
                    help="Name for the model")
solver.add_argument("--validate",dest="validate",
                    help="Validate solver results.")
log_group = parser.add_argument_group("log","For solver logging")
log_group.add_argument("--log",dest="log",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--console",dest="console",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--logfile",dest="logfile",action="store",
                       help="File to Log")

args = parser.parse_known_args()


class MIP(object):
    def __init__(self,instance=None):
        self.__instance=instance
        self.__model = None

    def build_model(self):
        pass

    def solve(self):
        return tuple([None,None,None])
    
    def lpsolve(self):
        return tuple([None,None,None])
        
class MIP1(MIP):
    def __init__(self,instance=None):
        super().__init__(instance=instance)
        
class MIP2(MIP):
    def __init__(self,instance=None):
        super().__init__(instance=instance)

class MIP3(MIP):
    def __init__(self,instance=None):
        super().__init__(instance=instance)

            
    
    def build_model(self):
        self.__build_model()

    def __build_model(self):
        if self.__model:
            return

        nres = self.__instance.nres
        nproc = self.__instance.nproc
        nmach = self.__instance.nmach
        
        R = self.__instance.R
        PMC = self.__instance.PMC
        bT = self.__instance.bT
        T = self.__instance.T
        C = self.__instance.C
        SC = self.__instance.SC
        MMC = self.__instance.MMC
        Wlc = self.__instance.Wlc
        Wbal = self.__instance.Wbal
        WPMC = self.__instance.WPMC
        WMMC = self.__instance.WMMC

        x_bar = self.__instance.map_assign()
        
        self.__model = Model()

        self.__model.Params.Quad = 1
        self.__model.Params.ScaleFlag = 0 
        self.__model.Params.NumericFocus = 3 

        self.__model.ModelSense = GRB.MINIMIZE

        self.__x = self.__model.addVars(nproc,nmach,
                                        vtype=GRB.BINARY,name="x")

        for p in range(nproc):
            for m in range(nmach):
                self.__x[p,m].start = x_bar[p,m]

        self.__z = self.__model.addVars(nproc,nmach,
                                        vtype=GRB.BINARY,name="z")
        
        
        self.__u = self.__model.addVars(nmach,nres,name="u",vtype=GRB.INTEGER,lb=0,ub=C)
        self.__u = self.__model.addVars(nmach,nres,name="ut",vtype=GRB.INTEGER,lb=0,ub=C)

        self.__t = self.__model.addVars(nmach,nmach,vtype=GRB.INTEGER,lb=0,name="t")

        self.__a = self.__model.addVars(nmach,nres,name="a",vtype=GRB.INTEGER,lb=0,ub=C)

        self.__d = self.__model.addVars(nmach,nres,name="d",vtype=GRB.INTEGER,lb=0,ub=C)
        
        self.__b = self.__model.addVars(nmach,nres,nres,vtype=GRB.INTEGER,name="b",lb=0)

        self.__o = self.__model.addVars(nserv,len(L),vtype=GRB.BINARY,name="o",lb=0)

        self.__g = self.__model.addVars(nserv,vtype=GRB.BINARY,name="s")

        self.__h = self.__model.addVars(nserv,len(N),vtype=GRB.BINARY,name="h")

        self.__lc = self.__model.addVars(nres,name="locadcost",vtype=GRB.INTEGER,lb=0,obj=Wlc)
        self.__bc= self.__model.addVars(nres,nres,name="balancecost",vtype=GRB.INTEGER,lb=0,obj=Wbal)
        self.__pmc = self.__model.addVar(name="pmc",vtype=GRB.INTEGER,lb=0,obj=WPMC)
        self.__smc = self.__model.addVar(name="smc",vtype=GRB.INTEGER,lb=0,obj=WSMC)
        self.__mmc = self.__model.addVar(name="mmc",vtype=GRB.INTEGER,lb=0,obj=WMMC)

        self.__model.update()

        self.__model.addConstrs((quicksum(self._x[p,m] for m in range(nmach)) == 1 for p in range(nproc)), name="process_assigned")

        self.__model.addConstrs((self.__u[m,r] == quicksum(self.__x[p,m]*R[p,r] for p in range(nproc)) for r in range(nres) for m in range(nmach)),name="util")
        self.__model.addConstrs((self.__ut[m,r] == quicksum(self.__z[p,m]*T[r]*R[p,r] for p in range(nproc)) for r in range(nres) for m in range(nmach)),name="transient")

        self.__model.addConstrs((self.__u[m,r]+self.__ut[m,r] <= C[m,r] for r in range(nres) for m in range(nmach)),name="cap")

        self.__model.addConstrs((self.__a[m,r] + self.__u[m,r] == C[m,r] for r in range(nres) for m in range(nmach)),name="avail")

        self.__model.addConstrs((self.__d[m,r] >= SC[m,r] - self.__u[m,r] for r in range(nres) for m in range(nmach)),name="load")

        self.__model.addConstrs((self.__z[p,m] >= x_bar[p,m] - self.__x[p,m] for p in range(nproc) for m in range(nmach)),name="z")

        self.__model.addConstrs((quicksum(self.__x[p,m] for p in s)<=1 for s in S for m in range(nmach)),name="conflict")

        self.__model.addConstrs((t[i,j]==quicksum(x_bar[p,i]*self.__x[p,j] for p in range(nproc)) for i in range(nmach)),name="t")

        self.__model.addConstrs((self.__o[s,l] >= self.__x[p,m] for p in s for m in l for s in S for l in L),name="o_lb")

        self.__model.addConstr((self.__o[s,l] <= quicksum((self.__x[p,m] for p in s for m in l)) for s in S for l in L),name="o_ub")

        self.__model.addConstr((quicksum(self.__o[s,l] for l in L) >= delta[s] for s in S),name="spread")

        

    def solve(self):
        if self.__model is None:
            raise Exception()

        nres = self.__instance.nres
        nproc = self.__instance.nproc
        nmach = self.__instance.nmach

        
        model = self.__model

        model.optimize()

        obj = model.ObjVal

        X = None
        alloc = None

        return tuple([obj,X,alloc])

    def lpsolve(self):
        if self.__model is None:
            raise Exception()

        nres = self.__instance.nres
        nproc = self.__instance.nproc
        nmach = self.__instance.nmach

        
        model = self.__model.relax()

        model.optimize()

        obj = model.ObjVal

        X = np.array([[self.__x[p,m].X for m in range(nmach)] for p in range(nproc)],dtype=np.float)
        alloc = np.array(nproc,dtype=np.int32)

        for p in range(nproc):
            for m in range(nmach):
                if self.__x[p,m] > 0.5:
                    alloc[p] = m
                    break
        
        return tuple([obj,X,alloc])
