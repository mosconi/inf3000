import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg3 import CG3
from gurobipy import *

def _cb4(model,where):
    if where == GRB.Callback.MIPSOL:
        nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
        if nodecnt > model._nproc:
            model.terminate()

class CG5(CG3):
        
    def __init__(self,instance,args):
        super().__init__(instance,args)
        self._cb = _cb4

    def lp2mip(self):
        self.__lp2mip()

    def build_lpmodel(self):
        super().build_lpmodel()
        self._lp.Params.Method = 2

    def build_column_model(self,machine):
        super().build_column_model(machine)
        self._mach[machine].model.Params.Method = 2
        self._mach[machine].model._nproc = self._instance.nproc


    def __lp2mip(self):
        self._mip = self._lp

        for v in self._mip.getVars():
            v.vtype = GRB.BINARY

        return 

        for m in range(self._instance.nmach):
            self._lbd[m,0].Start=1

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        S = self._instance.S
        N = self._instance.N
        iL = self._instance.iL
        iN = self._instance.iN
        WSMC = self._instance.WSMC

        x0 = self._instance.map_assign()
            
        self._h=self._mip.addVars(nserv,len(N),vtype=GRB.BINARY,
                                  lb=0,ub=1,name="h")

        self._o=self._mip.addVars(nserv,len(L),vtype=GRB.BINARY,
                                  lb=0,ub=1,name="o")
        
        self._mip.update()

        self._h_lb_constr=self._mip.addConstrs((0 - self._h[s,n] >= 0
                                                for s in sorted(S)
                                                for n in N
                                                for m in N[n]),
                                               name="h_lb_constr")
        self._h_ub_constr=self._mip.addConstrs((0 - self._h[s,n] >=0
                                                for s in sorted(S)
                                                for n in N),
                                               name="h_ub_constr")


        self._o_ub_constr=self._mip.addConstrs((-self._o[s,l] + 0 >=0
                                                for s in sorted(S)
                                                for l in L),
                                               name="o_ub_constr")
        self._o_lb_constr=self._mip.addConstrs((-self._o[s,l] + 0  <=0
                                                for s in sorted(S)
                                                for l in L
                                                for m in L[l]),
                                               name="o_lb_constr")

        self._dep_constr = self._mip.addConstrs((self._h[s,n] <= self._h[_s,n]
                                                 for n in N
                                                 for s in sdep
                                                 for _s in sdep[s]),
                                                name="dep")

        self._spread_constr = \
                              self._mip.addConstrs((self._o.sum(s,'*') >= delta[s]
                             for s in sorted(S)),
                                                   name="spread")

        
        self._mip.update()
        
        # for each column, complete the new constraints

        
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                procs = np.array(v._procs).flatten()
                z = procs - x0[:,m]
                z[z<0] = 0
                for s in sorted(S):
                     # é binário, por causa da restrição de conflito
                    serv = procs[S[s]].sum()
                    if serv >1:
                        print("s: %d, m: %d, sum: %d" % (s,m,serv))
                        print(v.VarName)
                        print(S[s])
                        for p in S[s]:
                            print(procs[p])
                    gserv = z[S[s]].sum()
                    if gserv >1:
                        print("s: %d, m: %d, gsum: %d" % (s,m,gserv))
                        print(S[s])
                    self._mip.chgCoeff( self._h_lb_constr[s,iN[m],m], v, serv)
                    self._mip.chgCoeff( self._h_ub_constr[s,iN[m]], v, serv)
                    self._mip.chgCoeff( self._o_lb_constr[s,iL[m],m], v, serv)
                    self._mip.chgCoeff( self._o_ub_constr[s,iL[m]], v, serv)

        self._mip.update()

    def presolve(self):
        self.__presolve()
        
    def __presolve(self):

        from collections import deque

        self._mip.Params.NumericFocus=1

        lbds = min([self._lbd.select(m,'*') for m in range(self._instance.nmach)],key=len)
        
        queue=deque()
        unqueue=deque()
        allvars = list()
        for m in range(self._instance.nmach):
            allvars.extend(self._lbd.select(m,'*'))
            
        for v in allvars:
            v._obj = v.obj
            v._star= False
            v.obj=0

        
        
        print("-"*80)
        print("  presolve %d" % (m))
        lbds[0]._star=1
        lbds[0].ub=0

        print("  presolve %d allvars" % (m))
        _allvars  = self._lbd.select(m,'*')
        queue.clear()
        
        v = None
        cond = True;
        while cond:
            print("  presolve %d optimize" % (m))
            self._mip.optimize()
            print("  presolve %d status: %d, len: %d" % (m,self._mip.status, len(queue)))
            if self._mip.status == 2:
                for _v in lbds:
                    if _v.X >.5:
                        v = _v
                for _v in allvars:
                    if _v.X >.5:
                        _v._star=True
            else:
                cond = False
                                
                        
            if v:
                print("  presolve %d len(queue): %d" % (m,len(queue)))
                v.ub = 0
                queue.append(v)
            else:
                cond = False

        for v in queue:
            v.ub=1

        lbds[0].ub=1
       
        counter = 0
        for v in allvars:
            v.obj = v._obj
            if v._star:
                counter +=1
            if not v._star :
                v.ub = 0

        print("counter: %d / %d" % (counter,len(allvars)))
        self._smc = self._mip.addVar(name="smc",vtype=GRB.INTEGER,
                                     lb=0,obj=self._instance.WSMC)
        self._g=self._mip.addVars(self._instance.nserv,vtype=GRB.BINARY,
                                  lb=0,ub=1,name="g")

        self._mip.update()
        S = self._instance.S
                
        self._g_constr = self._mip.addConstrs((-self._g[s] == 0 
                                              for s in sorted(S)),
                                              name="g")

        self._smc_constr = self._mip.addConstrs((self._smc >= self._g[s]
                                                for s in sorted(S)),
                                                name="smc")
        self._mip.update()
        x0 = self._instance.map_assign()
        
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                procs = np.array(v._procs).flatten()
                z = procs - x0[:,m]
                z[z<0] = 0
                for s in sorted(S):
                     # é binário, por causa da restrição de conflito
                    gserv = z[S[s]].sum()
                    if gserv >1:
                        print("s: %d, m: %d, gsum: %d" % (s,m,gserv))
                        print(S[s])
                    self._mip.chgCoeff( self._g_constr[s], v, gserv)

        self._mip.update()

        self._mip.Params.NumericFocus=3
