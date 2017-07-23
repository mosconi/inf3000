import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

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

        _tp = tuple(compcol.procs)
        for v in self._lbd.select(machine,"*"):
            if v._procs == _tp:
                return CGAdd(status=CGAddStatus.Exist,var=None)

        
        _col = Column()
        _col.addTerms(
            compcol.procs,
            [self._p_constr[p] for p in range(nproc)]
            )
        _col.addTerms(1,self._m_constr[machine])
        _col.addTerms(
            compcol.servs,
            [self._h_lb_constr[s,iN[machine],machine] for s in sorted(S)]
            )
        _col.addTerms(
            compcol.servs,
            [self._h_ub_constr[s,iN[machine]] for s in sorted(S)]
            )
        _col.addTerms(
            compcol.servs,
            [self._o_lb_constr[s,iL[machine],machine] for s in sorted(S)]
            )
        _col.addTerms(
            compcol.servs,
            [self._o_ub_constr[s,iL[machine]] for s in sorted(S)]
            )
        _col.addTerms(
            compcol.g,
            [self._g_constr[s] for s in sorted(S)]
            )

        _c = len(self._lbd.select(machine,'*'))
        var = self._lp.addVar(
            obj = compcol.obj,
            vtype=GRB.CONTINUOUS,
            ub=1,lb=0,
            column = _col,
            name = "lbd[%d,%d]" % (machine,len(self._lbd.select(machine,'*')))
        )
        self._lbd[machine,_c] = var
        var._procs = tuple(compcol.procs)
        #self._procs[machine][tuple(compcol.procs)] +=1
        
        return CGAdd(status=CGAddStatus.Added,var=var)
        
        
        
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
        self._lp.Params.Method = 1
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
            self._lbd[m,0]._procs = tuple(c)
            
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
                                            for s in sorted(S)
                                            for n in N
                                            for m in N[n]),
                                              name="h_lb_constr")
        self._h_ub_constr=self._lp.addConstrs((quicksum(self._lbd[m,0]*x0[p,m] for p in S[s] for m in N[n]) -self._h[s,n] >=0
                                            for s in sorted(S)
                                            for n in N),
                                              name="h_ub_constr")


        self._o_ub_constr=self._lp.addConstrs((-self._o[s,l] + quicksum(self._lbd[m,0]*x0[p,m] for p in S[s] for m in L[l]) >=0 
                                            for s in sorted(S)
                                            for l in L),
                                              name="o_ub_constr")
        self._o_lb_constr=self._lp.addConstrs((-self._o[s,l] + quicksum(self._lbd[m,0]*x0[p,m] for p in S[s])  <=0 
                                            for s in sorted(S)
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
                             for s in sorted(S)),
                                                  name="spread")

        self._g_constr = self._lp.addConstrs((-self._g[s] == 0 
                                              for s in sorted(S)),
                                             name="g")

        self._smc_constr = self._lp.addConstrs((self._smc >= self._g[s]
                                                for s in sorted(S)),
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
        T = self._instance.T
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
        gamma = np.ones(nserv,dtype=np.float64)   # duais de serviço migrado

        self._mach[machine].model = Model(name="mach_%d" % machine,env=self._env)

        self._mach[machine].model.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._mach[machine].model.Params.OutputFlag = 1
            self._mach[machine].model.Params.LogFile= self._args.logfile
            
        
        self._mach[machine].model.ModelSense = GRB.MINIMIZE

        self._mach[machine].model.Params.ScaleFlag = 0
        self._mach[machine].model.Params.Quad = 1
        self._mach[machine].model.Params.NumericFocus = 3

        self._mach[machine]._obj = self._mach[machine].model.addVar(vtype=GRB.INTEGER,name="obj",obj=0)

        self._mach[machine]._pixp = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="pixp",obj=0)
        self._mach[machine]._hsigma = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="hsigma",obj=0)
        self._mach[machine]._ggamma = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="ggamma",obj=0)


        
        self._mach[machine]._x = self._mach[machine].model.addVars(nproc,
                                                                   obj=pi,
                                                                   vtype=GRB.BINARY,name="x")

        # for p in range(nproc):
        #     self._mach[machine]._x[p].Start = x0[p]

        self._mach[machine]._z = self._mach[machine].model.addVars(nproc,ub=x0,
                                                                   vtype=GRB.BINARY,name="z")
        # for p in range(nproc):
        #     self._mach[machine]._z[p].Start = 0

        self._mach[machine]._g = self._mach[machine].model.addVars(nserv,obj=gamma,
                                                                   vtype=GRB.BINARY,name="g")  # serviço migrado
        # for s in sorted(S):
        #     self._mach[machine]._g[s].Start = 0

        self._mach[machine]._u = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="u")

        self._mach[machine]._ut = self._mach[machine].model.addVars(nres,
                                                                    vtype=GRB.INTEGER,lb=0,ub=C*T,name="ut")

        self._mach[machine]._d = self._mach[machine].model.addVars(nres,obj=Wlc,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="d")

        self._mach[machine]._a = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="a")

        self._mach[machine]._b = self._mach[machine].model.addVars(nres,nres,obj=Wbal,
                                                                   vtype=GRB.INTEGER,lb=0,name="b")

        # serve tanto para h[s,n], h[s,n,m] e o[s,l] o[s,l,m]
        self._mach[machine]._h = self._mach[machine].model.addVars(nserv,obj=sigma,
                                                                   vtype=GRB.BINARY,name="h") # serviço presente
#        for s in sorted(S):
#            self._mach[machine]._h[s].Start = x0[S[s]].sum()


        self._mach[machine]._pmc = self._mach[machine].model.addVar(vtype=GRB.INTEGER,
                                                                    lb=0,obj=WPMC,
                                                name="pmc")

        self._mach[machine]._mmc = self._mach[machine].model.addVar(vtype=GRB.INTEGER,
                                                                    lb=0,obj=WMMC,
                                                name="mmc")


        self._mach[machine]._cte = self._mach[machine].model.addVar(vtype=GRB.BINARY,
                                                obj=mu,
                                                lb=1,
                                                ub=1,
                                                name="Constant")

        self._mach[machine].model.update()

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._z[p] + self._mach[machine]._x[p] == 1
                for p in range(nproc)  if x0[p]==1),
            name="z"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._g[s] == quicksum(x0[p]*self._mach[machine]._z[p] for p in S[s])
                for s in sorted(S)),
            name="g"
        )
        # self._mach[machine].model.addConstrs(
        #     (
        #         self._mach[machine]._g[s] == self._mach[machine]._z.sum(S[s])
        #         for s in sorted(S)),
        #     name="g"
        # )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] == quicksum(self._mach[machine]._x[p]*R[p,r] for p in range(nproc))
                for r in range(nres)),
            name="util"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._ut[r] == quicksum(R[p,r]*self._mach[machine]._z[p] for p in range(nproc) if x0[p] == 1)*T[r]
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
                self._mach[machine]._a[r] ==  C[r] - self._mach[machine]._u[r]
                for r in range(nres)),
            name="avail"
        )
            
        self._mach[machine].model.addConstrs(
            (
                 self._mach[machine]._u[r] -  self._mach[machine]._d[r] <= SC[r]
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
                self._mach[machine]._x.sum(S[s])  <= 1
                for s in sorted(S)),
            name="conflict"
        )
        
        # self._mach[machine].model.addConstrs(
        #     (
        #         self._mach[machine]._h[s] >= self._mach[machine]._x[p] 
        #         for s in sorted(S)
        #         for p in S[s]
        #     ),
        #     name="h_lb"
        # )
        # self._mach[machine].model.addConstrs(
        #     (
        #         self._mach[machine]._h[s] <= self._mach[machine]._x.sum(S[s])
        #         for s in sorted(S)
        #     ),
        #     name="h_ub"
        # )
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._h[s] == self._mach[machine]._x.sum(S[s])
                for s in sorted(S)
            ),
            name="h"
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
            -self._mach[machine]._obj + 
            WMMC*self._mach[machine]._mmc +
            WPMC*self._mach[machine]._pmc +
            quicksum(Wlc[r]*self._mach[machine]._d[r] 
                     for r in range(nres))+
            quicksum(Wbal[r1,r2]*self._mach[machine]._b[r1,r2] 
                     for r1 in range(nres)
                     for r2 in range(nres)
                     )==0
            ,name="cost")

        self._mach[machine]._pixp_constr = self._mach[machine].model.addConstr(
            self._mach[machine]._pixp == 
            quicksum(pi[p]*self._mach[machine]._x[p] for p in range(nproc)) 
            ,name="pixp"
        )        
        self._mach[machine]._hsigma_constr = self._mach[machine].model.addConstr(
            self._mach[machine]._hsigma ==
            quicksum(sigma[s]*self._mach[machine]._h[s] for s in sorted(S)) 
            ,name="hsigma"
        )        
        self._mach[machine]._ggamma_constr = self._mach[machine].model.addConstr(
            self._mach[machine]._ggamma ==
            quicksum(gamma[s]*self._mach[machine]._g[s] for s in sorted(S)) 
            ,name="ggamma"
        )        

        self._mach[machine].model.update()

    def column_compute(self,machine):
        return self.__column_compute(machine)
    def column_prepare(self,machine, pi, mu, sigma, gamma):
        self.__column_prepare(machine, pi, mu, sigma, gamma)
    def __column_prepare(self,machine, pi, mu, sigma, gamma):
        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S

        x0 = self._instance.mach_map_assign(machine)

        self._mach[machine].model.reset()

        self._mach[machine]._cte.Obj = - mu
        
        for p in range(nproc):
            self._mach[machine]._x[p].Obj = - pi[p]
            if self._args.validate:
                pass
                # self._mach[machine].model.chgCoeff(
                #      self._mach[machine]._pixp_constr,
                #      self._mach[machine]._x[p],
                #      -pi[p]
                #      )


        for s in range(nserv):
            self._mach[machine]._h[s].Obj = - sigma[s]
            self._mach[machine]._g[s].Obj = - gamma[s]

            if self._args.validate:
                pass
                # self._mach[machine].model.chgCoeff(
                #     self._mach[machine]._ggamma_constr,
                #     self._mach[machine]._g[s],
                #     -gamma[s]
                # )
                # self._mach[machine].model.chgCoeff(
                #     self._mach[machine]._hsigma_constr,
                #     self._mach[machine]._h[s],
                #     -sigma[s]
                # )

        self._mach[machine].model.update()

        for p in range(nproc):
            self._mach[machine]._x[p].Start = x0[p]

    def column_write(self,machine,k):
        self.__column_write(machine,k)
        
    def __column_write(self,machine,k):
        if self._args.dump:
            self._mach[machine].model.write("%s_%d_%d.lp"%(self._args.run_name,machine, k))
            self._mach[machine].model.write("%s_%d_%d.mps"%(self._args.run_name,machine, k))
            self._mach[machine].model.write("%s_%d_%d.prm"%(self._args.run_name,machine, k))

    def column_writesol(self,machine,k):
        if self._args.dump:
            self.__column_writesol(machine,k)

    def __column_writesol(self,machine,k):
        if self._args.dump:
            self._mach[machine].model.write("%s_%d_%d.sol"%(self._args.run_name,machine, k))

        
    def __column_compute(self,machine):
            
        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S

        x0 = self._instance.mach_map_assign(machine)

        self._mach[machine].model.optimize()
        status = self._mach[machine].model.Status
        if status == GRB.UNBOUNDED:
            raise Exception("UNBOUNDED")

        elif status == GRB.INFEASIBLE:
            self._mach[machine].model.computeIIS()
            self._mach[machine].model.write("%s_%d.ilp"%(self._args.run_name,machine))
            raise Exception("INFEASIBLE")

        elif status == GRB.INF_OR_UNBD:
            self._mach[machine].model.computeIIS()
            self._mach[machine].model.write("%s_%d.ilp"%(self._args.run_name,machine))
            model = self._mach[machine].model.copy()
            __inf_rc = self.__inf(model)
            
            raise Exception("INF OR UNBD")
           
        q = np.array(
            [1* (self._mach[machine]._x[p].X > .5) for p in range(nproc)], dtype=np.int32
        )
        z = np.array(
            [1* (self._mach[machine]._z[p].X > .5) for p in range(nproc)], dtype=np.int32
        )

        serv = np.array(
            [1* (self._mach[machine]._h[s].X > .5) for s in sorted(S)], dtype=np.int32
        )
        
        g = np.array(
             [1* (self._mach[machine]._g[s].X > .5) for s in sorted(S)], dtype=np.int32
        )

        _d = None
        _u = None
        _ut = None
        _a = None
        _b = None
        _hsigma = None
        _ggamma = None
        _pixp = None
        _pmc = None
        _mmc = None
        
        if self._args.validate:
            _d = np.array([self._mach[machine]._d[r].X for r in range(nres)])
            _u = np.array([round(self._mach[machine]._u[r].X) for r in range(nres)],dtype=np.int32)
            _ut = np.array([round(self._mach[machine]._ut[r].X) for r in range(nres)],dtype=np.int32)
            _a = np.array([round(self._mach[machine]._a[r].X) for r in range(nres)],dtype=np.int32)
            _b = np.array([[round(self._mach[machine]._b[r1,r2].X) for r2 in range(nres)] for r1 in range(nres)],dtype=np.int32)
            _hsigma = self._mach[machine]._hsigma.X
            _ggamma = self._mach[machine]._ggamma.X
            _pixp = self._mach[machine]._pixp.X
            _pmc =self._mach[machine]._pmc.X
            _mmc = self._mach[machine]._mmc.X
            
        return CGColumn(rc=self._mach[machine].model.objVal,
                        rtime = self._mach[machine].model.Runtime,
                        obj=round(self._mach[machine]._obj.X),
                        procs=q,
                        servs=serv,
                        g= g,
                        z=z,
                        hsigma = _hsigma,
                        ggamma = _ggamma,
                        pixp = _pixp,
                        u = _u, 
                        ut = _ut,
                        a= _a,
                        d = _d,
                        b = _b,
                        pmc = _pmc,
                        mmc = _mmc
                        )
                      

    def solve_relax(self):
        return self.__solve_relax()

    def lpwrite(self,k):
        self._lp.write("%s_master_%d.lp"%(self._args.run_name, k))
        self._lp.write("%s_master_%d.mps"%(self._args.run_name, k))
        self._lp.write("%s_master_%d.prm"%(self._args.run_name, k))        

    def lpwritesol(self,k):
        self._lp.write("%s_master_%d.lp"%(self._args.run_name, k))

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
            self._g_constr[s].Pi for s in sorted(S)
        ],dtype=np.float64)

        return RelaxSolution(obj = _obj,
                             rtime = self._lp.RunTime,
                             pi = _pi,
                             mu = _mu,
                             gamma = _gamma,
                             omikron_lb = _omikron_lb,
                             omikron_ub = _omikron_ub,                         
                             eta_lb = _eta_lb,
                             eta_ub = _eta_ub
        )
    
                        
        


    def solve_lp(self):
        if not self._lp:
            raise Exception("LP model not defined")
        pass


    def column_validate(self, colres, machine,pi, mu, sigma, gamma):
        return self.__validate_column(colres, machine,pi, mu, sigma, gamma)

    def __validate_column(self, colres, machine,pi, mu, sigma, gamma):

        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        R = self._instance.R
        T = self._instance.T
        C = self._instance.C[machine]
        SC = self._instance.SC[machine]
        S = self._instance.S
        SM = self._instance.SM
        bT = self._instance.bT
        PMC = self._instance.PMC
        MMC = self._instance.MMC[self._instance.assign(),machine]
        Wlc = self._instance.Wlc
        Wbal = self._instance.Wbal
        WPMC = self._instance.WPMC
        WMMC = self._instance.WMMC

        x0 = self._instance.mach_map_assign(machine).reshape(1,nproc)

        _x = colres.procs.reshape(1,nproc)

        for s in sorted(S):
            if colres.procs[S[s]].sum() >1:
                print(" conflict failed!")
                return CGValidate(status=CGValidateStatus.Invalid)

        _h = _x.dot(SM)

        if np.any(_h!=colres.servs):
            print(" services doesn`t match")
            print(_h)
            print(colres.servs)
            print(_h - colres.servs)
            print(_x)
            for s in sorted(S):
                print(s,_h[0,s],colres.servs[s],colres.procs[S[s]].sum())
            return CGValidate(status=CGValidateStatus.Invalid)
            

        _z = x0 - _x
        _z[_z < 0] = 0

        if np.any(_z!=colres.z):
            print(" migration doesn`t match")
            print(_z)
            print(colres.z)
            print("x0")
            print(x0)
            print("x0 - x ")
            print(x0 -_x )
            print(_z - colres.z)
            return CGValidate(status=CGValidateStatus.Invalid)

        
        _g = _z.dot(SM)

        if colres.g is not None and np.any(_g!=colres.g):
            print(" migrate services doesn`t match")
            print(_g)
            print(colres.g)
            return CGValidate(status=CGValidateStatus.Invalid)
            
        
        
        _u = (_x.dot(R)).sum(axis=0)
        
        if np.any(abs(_u - colres.u) > self._args.epslon):
            print("u[r] mismatch")
            print(" u (calc):")
            print(_u)
            print(" u (grb):")
            print(colres.u)
            print(" u (delta):")
            print(_u - colres.u)
            if np.any(abs((_u - colres.u)*1.0/_u) > self._args.tol) and \
               not self._args.accept:
                return CGValidate(status=CGValidateStatus.CalcMismatch)
            else:
                print("Error too low, ignoring")
            
        _ut = (T*_z.dot(R)).sum(axis=0)
        
        if np.any(abs(_ut - colres.ut) > self._args.epslon):
            print("ut[r] mismatch")
            print(_ut)
            print(colres.ut)

        if any(_u + _ut >C):
            print("Over capacity!")
            return CGValidate(status=CGValidateStatus.Invalid)

        _a = C - _u
        if np.any(abs(_a - colres.a) > self._args.epslon):
            print("a[r] mismatch")
            print(" a (calc):")
            print(_a)
            print(" a (grb):")            
            print(colres.a)
            print(" a (delta):")            
            print(_a - colres.a)
            if np.any(abs((_u - colres.u)*1.0/_u) > self._args.tol):
                print("Error big ( err > tol)")
                if not self._args.accept:
                    return CGValidate(status=CGValidateStatus.CalcMismatch)
            else:
                print("Error too low, ignoring")
        
        _d = _u - SC
        _d[_d < 0 ] = 0

        if np.any(abs(_d - colres.d) > self._args.epslon):
            print("d[r] mismatch")
            print(" u (calc):")
            print(_u)
            print(" u (grb):")
            print(colres.u)
            print(" SC:")
            print(SC)
            print(" d (calc):")
            print(_d)
            print(" d (grb):")
            print(colres.d)
            print(" d (delta):")
            print(_d - colres.d)
            print(" WLc")
            print(Wlc)
            if np.any(abs((_u - colres.u)*1.0/_u) > self._args.tol):
                print("Error big ( err > tol)")
                if not self._args.accept:
                    return CGValidate(status=CGValidateStatus.CalcMismatch)
            else:
                print("Error too low, ignoring")
        
        _b = np.empty((nres,nres),dtype=np.int32)

        for r1 in range(nres):
            for r2 in range(nres):
                _b[r1,r2] = bT[r1,r2]*_a[r1] - _a[r2]

        _b[_b<0] = 0
        
        _pmc = (_x*PMC).reshape(nproc)
        _pmc[x0.reshape(nproc) == 1] = 0

        _mmc = (_x*MMC).reshape(nproc)
        _mmc[x0.reshape(nproc) == 1] = 0

        _obj = (_d*Wlc).sum() + (_b*Wbal).sum() +WPMC*(_pmc.sum()) + WMMC*(_mmc.sum())
        if abs(_obj - colres.obj) > self._args.epslon:        
            print("  OBJ")
            print("  obj (calc):")
            print(_obj)
            print("  obj (grb):")
            print(colres.obj)
            print("  obj (delta)")
            print(_obj - colres.obj)
            print("  max W")
            print( max([np.amax(Wbal),
                   np.amax(Wlc),
                   np.amax(WPMC),
                   np.amax(WMMC)]) )
            if  (_obj - colres.obj)*1.0/_obj > max([np.amax(Wbal),
                   np.amax(Wlc),
                   np.amax(WPMC),
                   np.amax(WMMC)]) * self._args.tol:
                print("Error big ( err > tol)")
                if not self._args.accept:
                    return CGValidate(status=CGValidateStatus.CalcMismatch)
            else:
                print("Error too low, ignoring")


        _pixp = _x.dot(pi).sum()
        if abs( _pixp - colres.pixp) > self._args.epslon and False:
            print("  pi*x")
            print("  pi*x (calc):")            
            print(_pixp)
            print("  pi*x (grb):")            
            print(colres.pixp)
            print("  pi*x (delta):")            
            print(_pixp - colres.pixp)
            if not self._args.accept:
                return CGValidate(status=CGValidateStatus.CalcMismatch)
        
        _hsigma = _h.dot(sigma).sum()
        if abs( _hsigma - colres.hsigma) > self._args.epslon and False:
            print("  h*sigma")
            print("  h*sigma (calc):")
            print(_hsigma)
            print("  h*sigma (calc):")
            print(colres.hsigma)
            print("  h*sigma (delta):")
            print(_hsigma - colres.hsigma)
            if not self._args.accept:
                return CGValidate(status=CGValidateStatus.CalcMismatch)

        _ggamma = _g.dot(gamma).sum()
        if abs( _ggamma - colres.ggamma) > self._args.epslon and False:
            print("  g*gamma")
            print("  g*gamma (calc):")            
            print(_ggamma)
            print("  g*gamma (grb):")
            print(colres.ggamma)
            print("  g*gamma (delta):")
            print(_ggamma - colres.ggamma)            
            if not self._args.accept:
                return CGValidate(status=CGValidateStatus.CalcMismatch)

        _rc = _obj - _pixp  - _hsigma  - _ggamma - mu
        if abs(_rc - colres.rc) > self._args.epslon:
            print("RC")
            print(" rc (calc):")
            print(_rc)
            print(" rc (grb):")
            print(colres.rc)
            print(" rc (delta):")
            print(_rc - colres.rc)
            print(" mu:")
            print(mu)
            print("  max W:")
            print( max([np.amax(Wbal),
                   np.amax(Wlc),
                   np.amax(WPMC),
                   np.amax(WMMC)]) )
            if  (_rc - colres.rc) > max([np.amax(Wbal),
                   np.amax(Wlc),
                   np.amax(WPMC),
                   np.amax(WMMC)]) * self._args.tol:
                print("Error big ( err > tol)")
                if not self._args.accept:
                    return CGValidate(status=CGValidateStatus.CalcMismatch)
            else:
                print("Error too low, ignoring")

        return CGValidate(status=CGValidateStatus.Valid)


    def lp2mip(self):
        self.__lp2mip()

    def __lp2mip(self):
        self._mip = self._lp

        for v in self._mip.getVars():
            v.vtype = GRB.INTEGER

        for m in range(self._instance.nmach):
            self._lbd[m,0].Start=1
            

    def __inf(self,model):

        orignumvars = model.NumVars

        model.feasRelaxS(0, False, True, False)

        model.optimize()
        status = model.status

        if status in (GRB.Status.INF_OR_UNBD, GRB.Status.INFEASIBLE, GRB.Status.UNBOUNDED):
            print('The relaxed model cannot be solved because it is infeasible or unbounded')
            return False
        if status != GRB.Status.OPTIMAL:
            print('Optimization was stopped with status %d' % status)
            return False

        print('\nSlack values:')

        model.write("Infeasible.lp")
        allvars = model.getVars()
        slacks = model.getVars()[orignumvars:]
        nslacks = model.getVars()[:orignumvars]
        constrs = model.getConstrs()
        v_constrs = []
        for sv in slacks:
            if sv.X > 1e-6:
                print('%s = %g' % (sv.VarName, sv.X))
                for c in constrs:
                    if model.getCoeff(c,sv) > 1e-6:
                        v_constrs.append(c)

        for c in v_constrs:
            print(c)
            for v in allvars:
                if model.getCoeff(c,v) > 1e-6:
                    print('%s = %g' % (v.VarName, v.X))
        # conflict[4]: x[8] + x[23] + x[64] + x[99] - h[4] + ArtP_conflict[4] - ArtN_conflict[4] = 0

#        for varname in [ "x[8]", "x[23]", "x[64]", "x[99]","h[4]" ]:
#            v = model.getVarByName(varname)
#            print('%s = %g' % (v.VarName, v.X))
            

        return True

        
