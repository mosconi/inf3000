
from .mip import MIP,Solution
import numpy as np
from gurobipy import *

class MIP2(MIP):
    def __init__(self,instance,args):
        super().__init__(instance,args)
        if args.logfile:
            self._env=Env(args.logfile)
        else:
            self._env=Env("mip2.log")
            
    def build_model(self):
        self.__build_model()


    def __build_model(self):

        self._mip = Model(name=self._args.run_name, env=self._env)

        self._mip.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._mip.Params.OutputFlag = 1
            self._mip.Params.LogFile= self._args.logfile
        
        
        self._mip.ModelSense = GRB.MINIMIZE

        self._mip.Params.ScaleFlag = 0
        self._mip.Params.Quad = 1
        self._mip.Params.NumericFocus = 3

        
        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        R = self._instance.R
        C = self._instance.C
        SC = self._instance.SC
        S = self._instance.S
        bT = self._instance.bT
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        N = self._instance.N
        PMC = self._instance.PMC
        MMC = self._instance.MMC
        Wlc = self._instance.Wlc
        Wbal = self._instance.Wbal
        WPMC = self._instance.WPMC
        WSMC = self._instance.WSMC
        WMMC = self._instance.WMMC

        x0 = self._instance.map_assign()

        self._x = self._mip.addVars(nproc,nmach,vtype=GRB.BINARY,name="x")

        for p in range(nproc):
            for m in range(nmach):
                self._x[p,m].Start = x0[p,m]

        if self._args.verbose:
            print(" creating vars")
        self._z =self._mip.addVars(nproc,nmach,vtype=GRB.BINARY,name="z")

        self._t = self._mip.addVars(nmach,nmach,vtype=GRB.INTEGER,lb=0,ub=nproc,name="t")

        self._u = self._mip.addVars(nmach,nres,vtype=GRB.INTEGER,lb=0,ub=C,name="u")
        self._ut = self._mip.addVars(nmach,nres,vtype=GRB.INTEGER,lb=0,ub=C,name="ut")
        self._a = self._mip.addVars(nmach,nres,vtype=GRB.INTEGER,lb=0,ub=C,name="a")
        self._d = self._mip.addVars(nmach,nres,vtype=GRB.INTEGER,lb=0,ub=C,name="d")
        self._b = self._mip.addVars(nmach,nres,nres,vtype=GRB.INTEGER,lb=0,name="b")

        self._o = self._mip.addVars(nserv,len(L),vtype=GRB.BINARY,name="o")
        self._g = self._mip.addVars(nserv,vtype=GRB.BINARY,name="g")

        self._h = self._mip.addVars(nserv,len(N),vtype=GRB.BINARY,name="h")

        self._lc = self._mip.addVars(nres,name="lc",vtype=GRB.INTEGER,lb=0,obj=Wlc)
        self._bc = self._mip.addVars(nres,nres,name="bc",vtype=GRB.INTEGER,lb=0,obj=Wbal)
        self._pmc = self._mip.addVar(name="pmc",vtype=GRB.INTEGER,lb=0,obj=WPMC)
        self._smc = self._mip.addVar(name="smc",vtype=GRB.INTEGER,lb=0,obj=WSMC)
        self._mmc = self._mip.addVar(name="mmc",vtype=GRB.INTEGER,lb=0,obj=WMMC)
        
        self._mip.update()
        if self._args.verbose:
            print(" creating constraints")

        if self._args.verbose>2:
            print(" creating constraints process assign")

        self._mip.addConstrs((self._x.sum(p,'*') == 1
                             for p in range(nproc)),
                             name="p_assign")

        if self._args.verbose>2:
            print(" creating constraints utilization")

        self._mip.addConstrs((self._u[m,r] == quicksum(self._x[p,m]*R[p,r]
                                                       for p in range(nproc))
                              for r in range(nres) for m in range(nmach)),
                             name="util")

        if self._args.verbose>2:
            print(" creating constraints transient utilization")

        self._mip.addConstrs((self._ut[m,r] == quicksum(self._z[p,m]*R[p,r]
                                                        for p in range(nproc))
                              for r in range(nres) for m in range(nmach)),
                             name="transient")
        
        if self._args.verbose>2:
            print(" creating constraints capacity")
        self._mip.addConstrs((self._u[m,r]+ self._ut[m,r] <= C[m,r]
                              for r in range(nres) for m in range(nmach)),
                             name="cap")

        if self._args.verbose>2:
            print(" creating constraints availibility")
        self._mip.addConstrs((self._a[m,r] + self._u[m,r] == C[m,r]
                              for r in range(nres) for m in range(nmach)),
                             name="avail")
        
        if self._args.verbose>2:
            print(" creating constraints load")
        self._mip.addConstrs((self._u[m,r] - self._d[m,r] <= SC[m,r]
                              for r in range(nres) for m in range(nmach)),
                             name="load")

        if self._args.verbose>2:
            print(" creating constraints z")
        self._mip.addConstrs((self._z[p,m] >= x0[p,m] - self._x[p,m]
                              for  p in range(nproc) for m in range(nmach)),
                             name="z")

        if self._args.verbose>2:
            print(" creating constraints conflict")
        self._mip.addConstrs((self._x.sum(S[s],m) <= 1
                              for s in S
                              for m in range(nmach)),
                             name="conflict")
        
        if self._args.verbose>2:
            print(" creating constraints t")
        self._mip.addConstrs((self._t[i,j] == quicksum(x0[p,i]*self._x[p,j]
                                                       for p in range(nproc))
                              for i in range(nmach)
                              for j in range(nmach)),
                             name="t")

        if self._args.verbose>2:
            print(" creating constraints o")
        self._mip.addConstrs((self._o[s,l] <= self._x.sum(S[s],L[l])
                              for s in S
                              for l in L),
                             name="o_ub")
        self._mip.addConstrs((self._o[s,l] >= self._x[p,m]
                              for s in S
                              for l in L
                              for p in S[s]
                              for m in L[l]),
                             name="o_lb")
        
        if self._args.verbose>2:
            print(" creating constraints spread")
        self._mip.addConstrs((self._o.sum(s,'*') >= delta[s]
                              for s in S),
                             name="spread")
        
        if self._args.verbose>2:
            print(" creating constraints g")
        self._mip.addConstrs((self._g[s] == self._z.sum(S[s],'*')
                              for s in S),
                             name="g")

        if self._args.verbose>2:
            print(" creating constraints h")
        self._mip.addConstrs((self._h[s,n] <= self._x.sum(S[s],N[n])
                              for s in S
                              for n in N),
                             name="h_ub")
        self._mip.addConstrs((self._h[s,n] >= self._x[p,m]
                              for s in S
                              for p in S[s]
                              for n in N
                              for m in N[n]),
                             name="h_lb") 

        if self._args.verbose>2:
            print(" creating constraints dependency")
        self._mip.addConstrs((self._h[s,n] <= self._h[_s,n]
                              for n in N
                              for s in sdep
                              for _s in sdep[s]),
                             name="dep")

        if self._args.verbose>2:
            print(" creating constraints b")
        self._mip.addConstrs((self._b[m,r1,r2] >= bT[r1,r2]*self._a[m,r1] - self._a[m,r2]
                              for m in range(nmach)
                              for r1 in range(nres)
                              for r2 in range(nres)
                              ),name="b")
        
        if self._args.verbose>2:
            print(" creating constraints lc")
        self._mip.addConstrs((self._lc[r] == self._d.sum('*',r)
                              for r in range(nres)
                              ),name="lc")

        if self._args.verbose>2:
            print(" creating constraints bc")
        self._mip.addConstrs((self._bc[r1,r2] == self._b.sum('*',r1,r2)
                              for r1 in range(nres)
                              for r2 in range(nres)                              
                              ),name="bc")

        if self._args.verbose>2:
            print(" creating constraints pmc")
        self._mip.addConstr(self._pmc == quicksum(self._z.sum(p,'*')*PMC[p]
                                                  for p in range(nproc)),
                            name="pmc")
        
        if self._args.verbose>2:
            print(" creating constraints smc")
        self._mip.addConstrs((self._smc >= self._g[s] for s in S), name="smc")

        if self._args.verbose>2:
            print(" creating constraints mmc")
        self._mip.addConstr(self._mmc == quicksum(MMC[i,j]*self._t[i,j]
                                                  for i in range(nmach)
                                                  for j in range(nmach)),
                            name="mmc")
        
        self._mip.update()


        
