import numpy as np
from .solution import RelaxSolution,CGColumn

from .cg import CG
from gurobipy import *


class CG2(CG):
    
    def __init__(self,instance,args):
        super().__init__(instance,args)

    def build_model(self):
        self.__build_model()

    def __build_model(self):
        self._mip 

    def lp_add_col(self,machine,compcol):

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S
        iL = self._instance.iL
        iN = self._instance.iN
        
        _col = Column()
        _col.addTerms(
            compcol.procs,
            [self._p_constr[p] for p in range(nproc)]
            )
        _col.addTerms(1,self._m_constr[machine])
        _col.addTerms(
            compcol.servs,
            [self._h_lb_constr[s,iN[machine],machine] for s in S]
            )
        _col.addTerms(
            compcol.servs,
            [self._h_ub_constr[s,iN[machine]] for s in S]
            )
        _col.addTerms(
            compcol.servs,
            [self._o_lb_constr[s,iL[machine],machine] for s in S]
            )
        _col.addTerms(
            compcol.servs,
            [self._o_ub_constr[s,iL[machine]] for s in S]
            )
        _col.addTerms(
            compcol.g,
            [self._g_constr[s] for s in S]
            )

        self._lbd[machine,len(self._lbd.select(machine,'*'))] = self._lp.addVar(
            obj = compcol.obj,
            vtype=GRB.CONTINUOUS,
            ub=1,lb=0,
            column = _col,
            name = "lbd[%d,%d]" % (machine,len(self._lbd.select(machine,'*')))
            )
        
        
    def build_lpmodel(self):
        if self._lp:
            return

        self._lp = Model(name="RM",env=self._env)

        self._lp.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._lp.Params.OutputFlag = 1
            self._lp.Params.LogFile= self._args.logfile
            
        
        self._lp.ModelSense = GRB.MINIMIZE

        self._lp.Params.ScaleFlag = 0
        self._lp.Params.Quad = 1
        self._lp.Params.NumericFocus = 3

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

        self._lbd=self._lp.addVars(nmach,1,vtype=GRB.CONTINUOUS,
                                   lb=0,ub=1,name="lbd")

        for m in range(nmach):
            c = x0[:,m]
            obj= self._instance.mach_objective(m,c)
            self._lbd[m,0].obj=obj
            
        self._g=self._lp.addVars(nserv,vtype=GRB.CONTINUOUS,
                                 lb=0,ub=1,name="g")


        self._h=self._lp.addVars(nserv,len(N),vtype=GRB.CONTINUOUS,
                                 lb=0,ub=1,name="h")

        self._o=self._lp.addVars(nserv,len(L),vtype=GRB.CONTINUOUS,
                                 lb=0,ub=1,name="o")

        self._smc = self._lp.addVar(name="smc",vtype=GRB.CONTINUOUS,
                                    lb=0,obj=WSMC)

        self._lp.update()

        self._p_constr=self._lp.addConstrs((quicksum(self._lbd[m,0]*x0[p,m] for m in range(nmach))==1 
                                            for p in range(nproc)),
                                           name="p_constr")
        self._m_constr=self._lp.addConstrs((self._lbd[m,0]==1 
                                            for m in range(nmach)),
                                           name="m_constr")

        self._h_lb_constr=self._lp.addConstrs((quicksum(self._lbd[m,0]*x0[p,m] for p in S[s]) - self._h[s,n] >= 0
                                            for s in S
                                            for n in N
                                            for m in N[n]),
                                              name="h_lb_constr")
        self._h_ub_constr=self._lp.addConstrs((quicksum(self._lbd[m,0]*x0[p,m] for p in S[s] for m in N[n]) -self._h[s,n] >=0
                                            for s in S
                                            for n in N),
                                              name="h_ub_constr")


        self._o_ub_constr=self._lp.addConstrs((-self._o[s,l] + quicksum(self._lbd[m,0]*x0[p,m] for p in S[s] for m in L[l]) >=0 
                                            for s in S
                                            for l in L),
                                              name="o_ub_constr")
        self._o_lb_constr=self._lp.addConstrs((-self._o[s,l] + quicksum(self._lbd[m,0]*x0[p,m] for p in S[s])  <=0 
                                            for s in S
                                            for l in L
                                            for m in L[l]),
                                              name="o_lb_constr")

        self._dep_constr = self._lp.addConstrs((self._h[s,n] <= self._h[_s,n]
                                                 for n in N
                                                 for s in sdep
                                                 for _s in sdep[s]),
                                               name="dep")

        self._spread_constr = \
                              self._lp.addConstrs((self._o.sum(s,'*') >= delta[s]
                             for s in S),
                                                  name="spread")

        self._g_constr = self._lp.addConstrs((-self._g[s] == 0 
                                              for s in S),
                                             name="g")

        self._smc_constr = self._lp.addConstrs((self._smc >= self._g[s]
                                                for s in S),
                                               name="smc")
        
        self._lp.update()
        
    def relax(self):
        pass

    def build_column_model(self,machine):
        self.__build_column_model(machine)

    def __build_column_model(self,machine):

        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        R = self._instance.R
        C = self._instance.C[machine]
        SC = self._instance.SC[machine]
        S = self._instance.S
        bT = self._instance.bT
        PMC = self._instance.PMC
        MMC = self._instance.MMC[self._instance.assign(),machine]
        Wlc = self._instance.Wlc
        Wbal = self._instance.Wbal
        WPMC = self._instance.WPMC
        WMMC = self._instance.WMMC

        x0 = self._instance.mach_map_assign(machine)

        pi = np.ones(nproc,dtype=np.float64)      # duais de processos
        mu = 1                                    # dual  de máquina
        sigma = np.ones(nserv,dtype=np.float64)   # soma das duais de serviços
        gamma = np.ones(nserv,dtype=np.float64) # duais de serviço migrado

        self._mach[machine].model = Model(env=self._env,name="mach_%d" % machine)

        self._mach[machine].model.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._mach[machine].model.Params.OutputFlag = 1
            self._mach[machine].model.Params.LogFile= self._args.logfile
            
        
        self._mach[machine].model.ModelSense = GRB.MINIMIZE

        self._mach[machine].model.Params.ScaleFlag = 0
        self._mach[machine].model.Params.Quad = 1
        self._mach[machine].model.Params.NumericFocus = 3

        self._mach[machine]._obj = self._mach[machine].model.addVar(vtype=GRB.INTEGER,name="obj",obj=1)
        self._mach[machine]._x = self._mach[machine].model.addVars(nproc,
                                                                   vtype=GRB.BINARY,name="x",obj=-pi)

        self._mach[machine]._z = self._mach[machine].model.addVars(nproc,
                                                                   vtype=GRB.BINARY,name="z")

        self._mach[machine]._g = self._mach[machine].model.addVars(nserv,
                                                                   vtype=GRB.BINARY,name="g",obj=-gamma)

        self._mach[machine]._u = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="u")

        self._mach[machine]._ut = self._mach[machine].model.addVars(nres,
                                                                    vtype=GRB.INTEGER,lb=0,ub=C,name="ut")

        self._mach[machine]._d = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="d")

        self._mach[machine]._a = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="a")

        self._mach[machine]._b = self._mach[machine].model.addVars(nres,nres,
                                                                   vtype=GRB.INTEGER,lb=0,name="b")

        self._mach[machine]._lc = self._mach[machine].model.addVars(nres,
                                                                    vtype=GRB.INTEGER,lb=0,name="lc")

        # serve tanto para h[s,n] e o[s,n]
        self._mach[machine]._s = self._mach[machine].model.addVars(nserv,
                                                                   vtype=GRB.BINARY,name="s",obj=-sigma)


        self._mach[machine]._pmc = self._mach[machine].model.addVar(vtype=GRB.INTEGER,
                                                lb=0,
                                                name="pmc")

        self._mach[machine]._mmc = self._mach[machine].model.addVar(vtype=GRB.INTEGER,
                                                lb=0,
                                                name="mmc")


        self._mach[machine]._cte = self._mach[machine].model.addVar(vtype=GRB.BINARY,
                                                obj=mu,
                                                lb=1,
                                                ub=1,
                                                name="Constant")

        self._mach[machine].model.update()

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._z[p] >= 
                x0[p] - self._mach[machine]._x[p]
                for p in range(nproc)),
            name="z"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._g[s] == self._mach[machine]._x.sum(S[s])
                for s in S),
            name="g"
        )



        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] == quicksum(R[p,r]*self._mach[machine]._x[p] for p in range(nproc))
                for r in range(nres)),
            name="util"
        )

        
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._ut[r] == quicksum(R[p,r]*self._mach[machine]._z[p] for p in range(nproc))
                for r in range(nres)),
            name="transient"
        )
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] + self._mach[machine]._ut[r] <= C[r]
                for r in range(nres)),
                name="capacity"
            )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._a[r] ==  C[r] - self._mach[machine]._ut[r]
                for r in range(nres)),
            name="avail"
        )
            
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._d[r] == self._mach[machine]._ut[r] -  SC[r]
                for r in range(nres)),
            name="load"
        )
        self._mach[machine].model.addConstrs(
                (
                    self._mach[machine]._b[r1,r2] >=  bT[r1,r2]*self._mach[machine]._a[r1] - self._mach[machine]._a[r2]
                for r1 in range(nres)
                for r2 in range(nres)),
                name="balance"
            )
            
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._x.sum(S[s]) <=1
                for s in S),
            name="conflict"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._x.sum(S[s]) <=1
                for s in S),
            name="conflict"
        )
        self._mach[machine].model.addConstr(
            self._mach[machine]._pmc == quicksum(PMC[p]*(1-x0[p])*self._mach[machine]._x[p]
                                                 for p in  range(nproc)),
            name="pmc"
        )
        self._mach[machine].model.addConstr(
            self._mach[machine]._mmc == quicksum(MMC[p]*(1-x0[p])*self._mach[machine]._x[p]
                                                 for p in  range(nproc)),
            name="mmc"
        )
        self._mach[machine].model.addConstr(
            self._mach[machine]._obj == 
            self._mach[machine]._mmc +
            self._mach[machine]._pmc +
            quicksum(Wlc[r]*self._mach[machine]._lc[r] 
                     for r in range(nres))+
            quicksum(Wbal[r1,r2]*self._mach[machine]._b[r1,r2] 
                     for r1 in range(nres)
                     for r2 in range(nres))
            ,name="roadef_obj")

        self._mach[machine].model.update()

    def compute_column(self,machine, pi, mu, sigma, gamma):
        return self.__compute_column(machine, pi, mu, sigma, gamma)

    def __compute_column(self,machine, pi, mu, sigma, gamma):
            
        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv

        for p in range(nproc):
            self._mach[machine]._x[p].Obj = -pi[p]

        self._mach[machine]._cte.Obj = - mu

        for s in range(nserv):
            self._mach[machine]._s[s].Obj = -sigma[s]
            self._mach[machine]._g[s].Obj = -gamma[s]

        self._mach[machine].model.update()

        self._mach[machine].model.reset()

        self._mach[machine].model.optimize()

        q = np.array(
            [1* (self._mach[machine]._x[p].X > .5) for p in range(nproc)], dtype=np.int32
        )

        serv = np.array(
            [1* (self._mach[machine]._s[s].X > .5) for s in range(nserv)], dtype=np.int32
        )
        g = np.array(
            [1* (self._mach[machine]._g[s].X > .5) for s in range(nserv)], dtype=np.int32
        )

        return CGColumn(rc=self._mach[machine].model.objVal,
                      obj=round(self._mach[machine]._obj.X),
                      procs=q,
                      servs=serv,
                      g=g)
                      

    def solve_relax(self):
        return self.__solve_relax()

    def __solve_relax(self):
        if not self._lp:
            raise Exception("LP model not defined")

        self._lp.optimize()

        s = self._lp.Status
        if s == GRB.Status.INF_OR_UNBD:
            print("inf or unbd")
        elif s == GRB.Status.INFEASIBLE:
            self.lpiis()
        if s == GRB.Status.INF_OR_UNBD or \
           s == GRB.Status.INFEASIBLE or \
           s == GRB.Status.UNBOUNDED or \
           False:
            print("sdx")
            return None

        _obj = self._lp.objVal

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S
        N = self._instance.N
        L = self._instance.L        
        iN = self._instance.iN
        iL = self._instance.iL        

        _pi = np.array([self._p_constr[p].Pi for p in range(nproc)]
                       , dtype=np.float64)
        _mu = np.array([self._m_constr[m].Pi for m in range(nmach)]
                       , dtype=np.float64)
        _eta_lb = np.array([
            [self._h_lb_constr[s,iN[m],m].Pi for m in range(nmach) ]
            for s in range(nserv) ]
                           , dtype=np.float64)
        _eta_ub = np.array([
            [self._h_ub_constr[s,n].Pi for n in N ]
            for s in range(nserv) ]
                           , dtype=np.float64)
        _omikron_lb = np.array([
            [self._o_lb_constr[s,iL[m],m].Pi for m in range(nmach) ]
            for s in range(nserv) ]
                           , dtype=np.float64)
        _omikron_ub = np.array([
            [self._o_ub_constr[s,l].Pi for l in L ]
            for s in range(nserv) ]
                           , dtype=np.float64)
        _gamma = np.array([
            self._g_constr[s].Pi for s in S
        ],dtype=np.float64)

        return RelaxSolution(obj = _obj,
                         pi = _pi,
                         mu = _mu,
                         gamma = _gamma,
                         omikron_lb = _omikron_lb,
                         omikron_ub = _omikron_ub,                         
                         eta_lb = _eta_lb,
                         eta_ub = _eta_ub
                         )
                         
                        
        

    def solve(self):
        pass

    def solve_lp(self):
        if not self._lp:
            raise Exception("LP model not defined")
        pass

