import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg import CG
from gurobipy import *


class CG3(CG):
    
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
        # _col.addTerms(
        #     compcol.servs,
        #     [self._h_lb_constr[s,iN[machine],machine] for s in S]
        #     )
        # _col.addTerms(
        #     compcol.servs,
        #     [self._h_ub_constr[s,iN[machine]] for s in S]
        #     )
        # _col.addTerms(
        #     compcol.servs,
        #     [self._o_lb_constr[s,iL[machine],machine] for s in S]
        #     )
        # _col.addTerms(
        #     compcol.servs,
        #     [self._o_ub_constr[s,iL[machine]] for s in S]
        #     )
        # _col.addTerms(
        #     compcol.g,
        #     [self._g_constr[s] for s in S]
        #     )

        _c = len(self._lbd.select(machine,'*'))
        var = self._lp.addVar(
            obj = compcol.obj,
            vtype=GRB.CONTINUOUS,
            ub=1,lb=0,
            column = _col,
            name = "lbd[%d,%d]" % (machine,len(self._lbd.select(machine,'*')))
        )
        self._lbd[machine,_c] = var
        var._procs = compcol.procs
        #self._procs[machine][tuple(compcol.procs)] +=1
        
        return CGAdd(status=CGAddStatus.Added,var =var)
        
        
        
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
            self._lbd[m,0]._procs = c 
            
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

        self._mach[machine]._pixp = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="pixp",obj=-1)
        self._mach[machine]._hsigma = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="hsigma",obj=-1)
        self._mach[machine]._ggamma = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="ggamma",obj=-1)


        
        self._mach[machine]._x = self._mach[machine].model.addVars(nproc,
                                                                   vtype=GRB.BINARY,name="x")

        for p in range(nproc):
            self._mach[machine]._x[p].Start = x0[p]

        self._mach[machine]._u = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="u")

        self._mach[machine]._ut = self._mach[machine].model.addVars(nres,
                                                                    vtype=GRB.INTEGER,lb=0,ub=C*T,name="ut")

        self._mach[machine]._d = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="d")

        self._mach[machine]._a = self._mach[machine].model.addVars(nres,
                                                                   vtype=GRB.INTEGER,lb=0,ub=C,name="a")

        self._mach[machine]._b = self._mach[machine].model.addVars(nres,nres,
                                                                   vtype=GRB.INTEGER,lb=0,name="b")


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

#        self._mach[machine].model.addConstrs(
#            (
#                -self._mach[machine]._g[s] == -self._mach[machine]._z.sum(S[s])
#                for s in S),
#            name="g"
#        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] == quicksum(R[p,r]*self._mach[machine]._x[p] for p in range(nproc))
                for r in range(nres)),
            name="util"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._ut[r] == T[r]*quicksum(R[p,r]*(1-self._mach[machine]._x[p]) for p in range(nproc) if x0[p]==1)
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
                 self._mach[machine]._d[r] >= self._mach[machine]._u[r] -  SC[r]
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
                for s in S),
            name="conflict"
        )
        

        self._mach[machine].model.addConstr(
            self._mach[machine]._pmc == quicksum(PMC[p]*self._mach[machine]._x[p]
                                                 for p in  range(nproc) if x0[p] == 0 ),
            name="pmc"
        )
        self._mach[machine].model.addConstr(
            self._mach[machine]._mmc == quicksum(MMC[p]*self._mach[machine]._x[p]
                                                 for p in  range(nproc) if x0[p] == 0),
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
            -self._mach[machine]._pixp +
            quicksum(pi[p]*self._mach[machine]._x[p] for p in range(nproc)) ==0
            ,name="pixp"
        )        
        self._mach[machine]._hsigma_constr = self._mach[machine].model.addConstr(
            -self._mach[machine]._hsigma +
            quicksum(sigma[s]*self._mach[machine]._x[p] for s in S for p in S[s]) == 0
            ,name="hsigma"
        )        
        self._mach[machine]._ggamma_constr = self._mach[machine].model.addConstr(
            -self._mach[machine]._ggamma +
            quicksum(gamma[s]*(1-self._mach[machine]._x[p]) for s in S for p in S[s] if x0[p]==1) == 0
            ,name="ggamma"
        )        

        self._mach[machine].model.update()

    def compute_column(self,machine, pi, mu, sigma, gamma):
        return self.__compute_column(machine, pi, mu, sigma, gamma)

    def __compute_column(self,machine, pi, mu, sigma, gamma):
            
        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S
        x0 = self._instance.mach_map_assign(machine)

        for p in range(nproc):
            self._mach[machine].model.chgCoeff(
                self._mach[machine]._pixp_constr,
                self._mach[machine]._x[p],
                pi[p]
                )

        self._mach[machine]._cte.Obj = - mu

        g_acc =0
        for s in S:
            for p in S[s]:
                if x0[p] ==1:
                    self._mach[machine].model.chgCoeff(
                        self._mach[machine]._ggamma_constr,
                        self._mach[machine]._x[p],
                        -gamma[s]
                    )
                    g_acc += gamma[s]
                self._mach[machine].model.chgCoeff(
                    self._mach[machine]._hsigma_constr,
                    self._mach[machine]._x[p],
                    sigma[s]
                )
             
        self._mach[machine]._ggamma_constr = g_acc

        self._mach[machine].model.update()

        self._mach[machine].model.reset()
        self._mach[machine].model.write("%s_%d.lp"%(self._args.run_name,machine))
        self._mach[machine].model.write("%s_%d.mps"%(self._args.run_name,machine))

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
        
        self._mach[machine].model.write("%s_%d.sol"%(self._args.run_name,machine))
           
        q = np.array(
            [1* (self._mach[machine]._x[p].X > .5) for p in range(nproc)], dtype=np.int32
        )

        z = np.zeros(nproc,dtype = np.int32)
        for p in range(nproc):
            if x0[p] == 1:
                z[p] = round(1-self._mach[machine]._x[p].X)
        # z = np.array(
        #     [1* (self._mach[machine]._z[p].X > .5) for p in range(nproc)], dtype=np.int32
        # )

        serv = np.zeros(nserv,dtype = np.int32)
        for s in S:
            serv[s] = sum([round(1-self._mach[machine]._x[p].X) for p in S[s]])
       
        # serv = np.array(
        #     [1* (self._mach[machine]._x[p].X > .5) for s in S], dtype=np.int32
        # )

        g = np.zeros(nserv,dtype = np.int32)
        for s in S:
            for p in S[s]:
                if x0[p] == 1:
                    g[s] = round(1-self._mach[machine]._x[p].X)
                    break
        # g = np.array(
        #     [1* (self._mach[machine]._g[s].X > .5) for s in S], dtype=np.int32
        # )
        _d = np.array([self._mach[machine]._d[r].X for r in range(nres)])
        return CGColumn(rc=self._mach[machine].model.objVal,
                        obj=round(self._mach[machine]._obj.X),
                        procs=q,
                        servs=serv,
                        g=g,
                        z=z,
                        hsigma = self._mach[machine]._hsigma.X,
                        ggamma = self._mach[machine]._ggamma.X,
                        pixp = self._mach[machine]._pixp.X,
                        u = np.array([round(self._mach[machine]._u[r].X) for r in range(nres)],dtype=np.int32),
                        ut = np.array([round(self._mach[machine]._ut[r].X) for r in range(nres)],dtype=np.int32),
                        
                        a = np.array([round(self._mach[machine]._a[r].X) for r in range(nres)],dtype=np.int32),
                        
                        d = _d,
                        b = np.array([[round(self._mach[machine]._b[r1,r2].X) for r2 in range(nres)] for r1 in range(nres)],dtype=np.int32),
                        pmc = self._mach[machine]._pmc.X,
                        mmc = self._mach[machine]._mmc.X
                        )
                      

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
                         
                        
        


    def solve_lp(self):
        if not self._lp:
            raise Exception("LP model not defined")
        pass


    def validate_column(self, colres, machine,pi, mu, sigma, gamma):
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

        for s in S:
            if colres.procs[S[s]].sum() >1:
                print("conflict!")
                return CGValidate(status=CGValidateStatus.Invalid)

        _h = _x.dot(SM)

        if np.any(_h!=colres.servs):
            print(" services doesn`t match")
            print(_h)
            print(colres.servs)
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

        if np.any(_g!=colres.g):
            print(" services doesn`t match")
            print(_g)
            print(colres.g)
            return CGValidate(status=CGValidateStatus.Invalid)
            
        
        
        _u = (_x.dot(R)).sum(axis=0)
        
        if np.any(abs(_u - colres.u) > self._args.epslon):
            print("u[r] mismatch")
            print(_u)
            print(colres.u)
            
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
            print(_a)
            print(colres.a)

        
        _d = _u - SC
        _d[_d < 0 ] = 0

        if np.any(abs(_d - colres.d) > self._args.epslon):
            print("d[r] mismatch")
            print(_u)
            print(colres.u)
            print(SC)
            print(_d)
            print(colres.d)
            print(Wlc)
        
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
        _pixp = _x.dot(pi).sum()
        _hsigma = _h.dot(sigma).sum()
        _ggamma = _g.dot(gamma).sum()
        _rc = _obj - _pixp  - _hsigma  - _ggamma - mu
        
        if abs(_rc - colres.rc) > self._args.epslon:
            print("RC")
            print(_rc)
            print(colres.rc)
            print(mu)
            if abs(_obj - colres.obj) > self._args.epslon:
                print("obj")
                print(_obj)
                print(colres.obj)
            if abs( _hsigma - colres.hsigma) > self._args.epslon:
                print("h*sigma")
                print(_hsigma)
                print(colres.hsigma)
            if abs( _pixp - colres.pixp) > self._args.epslon:
                print("pi*x")
                print(_pixp)
                print(colres.pixp)
                
            if abs( _ggamma - colres.ggamma) > self._args.epslon:
                print("g*gamma")
                print(_ggamma)
                print(colres.ggamma)
            
            return CGValidate(status=CGValidateStatus.CalcMismatch)

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

        model.feasRelaxS(0, False, False, True)

        model.optimize()
        status = model.status

        if status in (GRB.Status.INF_OR_UNBD, GRB.Status.INFEASIBLE, GRB.Status.UNBOUNDED):
            print('The relaxed model cannot be solved \
            because it is infeasible or unbounded')
            return False
        if status != GRB.Status.OPTIMAL:
            print('Optimization was stopped with status %d' % status)
            return False

        print('\nSlack values:')

        model.write("Infeasible.lp")
        slacks = model.getVars()[orignumvars:]
        for sv in slacks:
            if sv.X > 1e-6:
                print('%s = %g' % (sv.VarName, sv.X))
        # conflict[4]: x[8] + x[23] + x[64] + x[99] - h[4] + ArtP_conflict[4] - ArtN_conflict[4] = 0

#        for varname in [ "x[8]", "x[23]", "x[64]", "x[99]","h[4]" ]:
#            v = model.getVarByName(varname)
#            print('%s = %g' % (v.VarName, v.X))
            

        return True

        
