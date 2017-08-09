import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg1 import CG1
from gurobipy import *


class CG3(CG1):
        
    def __init__(self,instance,args):
        super().__init__(instance,args)

    def lp2mip(self):
        self.__lp2mip()

    def build_lpmodel(self):
        super().build_lpmodel()
        self._lp.Params.Method = 2

    def build_column_model(self,machine):
        super().build_column_model(machine)
        #self._mach[machine].model.Params.Method = 2


    def __lp2mip(self):
        self._mip = self._lp

        for v in self._mip.getVars():
            v.vtype = GRB.INTEGER

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
            
        self._g=self._mip.addVars(nserv,vtype=GRB.BINARY,
                                  lb=0,ub=1,name="g")


        self._h=self._mip.addVars(nserv,len(N),vtype=GRB.BINARY,
                                  lb=0,ub=1,name="h")

        self._o=self._mip.addVars(nserv,len(L),vtype=GRB.BINARY,
                                  lb=0,ub=1,name="o")

        self._smc = self._mip.addVar(name="smc",vtype=GRB.INTEGER,
                                     lb=0,obj=WSMC)
        
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

        self._g_constr = self._mip.addConstrs((-self._g[s] == 0 
                                              for s in sorted(S)),
                                              name="g")

        self._smc_constr = self._mip.addConstrs((self._smc >= self._g[s]
                                                for s in sorted(S)),
                                                name="smc")
        
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
                    gserv = z[S[s]].sum()
                    self._mip.chgCoeff( self._h_lb_constr[s,iN[m],m], v, serv)
                    self._mip.chgCoeff( self._h_ub_constr[s,iN[m]], v, serv)
                    self._mip.chgCoeff( self._o_lb_constr[s,iL[m],m], v, serv)
                    self._mip.chgCoeff( self._o_ub_constr[s,iL[m]], v, serv)
                    self._mip.chgCoeff( self._g_constr[s], v, gserv)

        self._mip.update()
                                  
                    
