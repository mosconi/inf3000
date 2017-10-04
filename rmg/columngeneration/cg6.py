import numpy as np
from .solution import CGSolution,RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from gurobipy import *

class CG6(object):
    _instance = None
    _args = None
    _model = None
    _env = None
    
    def __init__(self,instance,args):
        self._args = args
        self._instance = instance

        nproc = self._instance.nproc
        P = self._instance.P
        nmach = self._instance.nmach
        M = self._instance.M
        
        self._cb = lambda model, where: None
        self._mach = [lambda: None for m in M]
        print([m for m in M])
        self._lbd = [[] for p in P]

        self._ppcounts = np.zeros((nproc, nproc), dtype=np.int32)
        self._mpcounts = np.zeros((nmach, nproc), dtype=np.int32)



    def __del__(self):
        M = self._instance.M

        if self._model:
            del(self._model)
            
        del(self._env)

    def load(self,filename):
        self._env = Env(self._args.logfile)
        self._model = read(filename,self._env)

    def save(self,filename):
        self._model.write(filename)

    def build_model(self):
        self._env=Env(self._args.logfile)
        if self._model:
            return

        self._model = Model(name="RM",env=self._env)

        self._model.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._model.Params.OutputFlag = 1
            self._model.Params.LogFile= self._args.logfile
            
        
        self._model.ModelSense = GRB.MINIMIZE

        self._model.Params.ScaleFlag = 0
        #self._model.Params.Method = 1
        self._model.Params.Quad = 1
        self._model.Params.NumericFocus = 3

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        M = self._instance.M
        P = self._instance.P
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

        for m in M:
            mval = self._instance.mach_validate(m,x0[:,m])
            var = self._model.addVar(vtype=GRB.CONTINUOUS, obj= mval.obj,
                                                   lb=0,ub=1,name="lbd[%d,0]"%(m))
            var._procs = x0[:,m]
            self._lbd[m].append(var)
            
        
        self._model.update()

        self._p_constr=self._model.addConstrs((quicksum(self._lbd[m][0]*x0[p,m] for m in M)==1 
                                               for p in P),
                                              name="p_constr")
        self._m_constr=self._model.addConstrs((self._lbd[m][0]==1 
                                               for m in M),
                                              name="m_constr")

        
        self._model.update()
        self._z_int = np.inf

    
    def build_column_model(self,machine):

        nproc = self._instance.nproc
        nres = self._instance.nres
        nserv = self._instance.nserv
        R = self._instance.R
        SR = self._instance.SR
        T = self._instance.T
        C = self._instance.C[machine]
        SC = self._instance.SC[machine]
        S = self._instance.S
        P = self._instance.P
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

        self._mach[machine].model = Model(name="mach_%d" % machine,env=self._env)

        self._mach[machine].model.Params.LogToConsole = self._args.console

        if self._args.logfile:
            self._mach[machine].model.Params.OutputFlag = 1
            self._mach[machine].model.Params.LogFile= self._args.logfile
            
        
        self._mach[machine].model.ModelSense = GRB.MINIMIZE

        self._mach[machine].model.Params.ScaleFlag = 0
        self._mach[machine].model.Params.Quad = 1
        self._mach[machine].model.Params.NumericFocus = 3
        self._mach[machine].model.Params.MIPGap = 0.01

        self._mach[machine]._obj = self._mach[machine].model.addVar(vtype=GRB.INTEGER,name="obj",obj=0)

        if self._args.validate:
            self._mach[machine]._pixp = self._mach[machine].model.addVar(vtype=GRB.CONTINUOUS,name="pixp",obj=0)

        self._mach[machine]._x = self._mach[machine].model.addVars(nproc,
                                                                   obj=pi,
                                                                   vtype=GRB.BINARY,name="x")

        # for p in P:
        #     self._mach[machine]._x[p].Start = x0[p]

        self._mach[machine]._z = self._mach[machine].model.addVars(nproc,ub=x0,
                                                                   vtype=GRB.BINARY,name="z")
        # for p in P:
        #     self._mach[machine]._z[p].Start = 0

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
                for p in P  if x0[p]==1),
            name="z"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] == quicksum(self._mach[machine]._x[p]*R[p,r] for p in P)
                for r in SR ),
            name="util"
        )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._ut[r] == quicksum(R[p,r]*self._mach[machine]._z[p] for p in P if x0[p] == 1)*T[r]
                for r in SR ),
            name="transient"
        )
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._u[r] + self._mach[machine]._ut[r] <= C[r]
                for r in SR),
                name="capacity"
            )

        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._a[r] ==  C[r] - self._mach[machine]._u[r]
                for r in SR),
            name="avail"
        )
            
        self._mach[machine].model.addConstrs(
            (
                 self._mach[machine]._u[r] -  self._mach[machine]._d[r] <= SC[r]
                for r in SR),
            name="load"
        )
        self._mach[machine].model.addConstrs(
                (
                    self._mach[machine]._b[r1,r2] >=  bT[r1,r2]*self._mach[machine]._a[r1] - self._mach[machine]._a[r2]
                    for r1 in SR
                    for r2 in SR),
                name="balance"
            )
            
        self._mach[machine].model.addConstrs(
            (
                self._mach[machine]._x.sum(S[s])  <= 1
                for s in sorted(S)),
            name="conflict"
        )
        
        self._mach[machine].model.addConstr(
            self._mach[machine]._pmc == quicksum(PMC[p]*(1-x0[p])*self._mach[machine]._x[p]
                                                 for p in  P),
            name="pmc"
        )
        self._mach[machine].model.addConstr(
            self._mach[machine]._mmc == quicksum(MMC[p]*(1-x0[p])*self._mach[machine]._x[p]
                                                 for p in  P),
            name="mmc"
        )
        self._mach[machine].model.addConstr(
            -self._mach[machine]._obj + 
            WMMC*self._mach[machine]._mmc +
            WPMC*self._mach[machine]._pmc +
            quicksum(Wlc[r]*self._mach[machine]._d[r] 
                     for r in SR)+
            quicksum(Wbal[r1,r2]*self._mach[machine]._b[r1,r2] 
                     for r1 in SR
                     for r2 in SR
                     )==0
            ,name="cost")

        if self._args.validate:
            self._mach[machine]._pixp_constr = self._mach[machine].model.addConstr(
                self._mach[machine]._pixp == 
                quicksum(pi[p]*self._mach[machine]._x[p] for p in P) 
                ,name="pixp"
            )        

        self._mach[machine].model.update()
        

    def write(self):
        if not self._model:
            raise Exception("Model not defined")

        self._model.write("%s.lp" %self._args.run_name)
        self._model.write("%s.mps" %self._args.run_name)



    def iis(self):
        if not self._model:
            raise Exception("Model not defined")
        
        self._model.computeIIS()
        self._model.write("%s.ilp" %self._args.run_name)


    def final_solve(self):
        if not self._model:
            raise Exception("Model not defined")

        self._model.optimize()
        status = self._model.Status
        if status == GRB.UNBOUNDED:
            raise Exception("UNBOUNDED")

        elif status == GRB.INFEASIBLE:
            self.iis()
            raise Exception("INFEASIBLE")

        elif status == GRB.INF_OR_UNBD:
            self.iis()
            raise Exception("INF OR UNBD")

                    

        _obj = self._model.ObjVal

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        P = self._instance.P
        M = self._instance.M

        x = np.zeros((nproc,nmach),dtype=np.int32)
        assign = np.empty(nproc,dtype=np.int32)
        for m in M:
            v = [ v for v in self._lbd[m]  if v.X>0 ][0]
            p = np.array([
                int(round(self._model.getCoeff(self._p_constr[p0], v)))
                for p0 in P
            ])
            x[:,m] = p
            assign[p==1] = m

        return CGSolution(obj = _obj, X=x, assign = assign)
        
    
    def solve(self):
        if not self._model:
            raise Exception("Model not defined")

        self._model.optimize()

        s = self._model.Status
        if s == GRB.Status.INF_OR_UNBD:
            self.lpiis()
            print("inf or unbd")
        elif s == GRB.Status.INFEASIBLE:
            self.lpiis()
            print("Infeasible ")
        elif s == GRB.Status.UNBOUNDED:
            print("unbounded")
        if s == GRB.Status.INF_OR_UNBD or \
           s == GRB.Status.INFEASIBLE or \
           s == GRB.Status.UNBOUNDED or \
           False:
            return None

        _obj = self._model.objVal

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        P = self._instance.P
        M = self._instance.M
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S
        N = self._instance.N
        L = self._instance.L        
        iN = self._instance.iN
        iL = self._instance.iL        

        _pi = np.array([self._p_constr[p].Pi for p in P]
                       , dtype=np.float64)
        _mu = np.array([self._m_constr[m].Pi for m in M]
                       , dtype=np.float64)
        _eta_lb = np.array([
            [0  for m in M ]
            for s in S ]
                           , dtype=np.float64)
        
        _eta_ub = np.array([
            [0 for n in N ]
            for s in S ]
                           , dtype=np.float64)
        _omikron_lb = np.array([
            [0 for m in M ]
            for s in S ]
                           , dtype=np.float64)
        _omikron_ub = np.array([
            [0 for l in L ]
            for s in S ]
                           , dtype=np.float64)
        _gamma = np.array([
            0 for s in sorted(S)
        ],dtype=np.float64)

        if abs(_obj - int(round(_obj)))< self._args.tol:
            _allint = all([[abs(round(v.x) - v.x ) < self._args.tol for v in self._lbd[m]] for m in M ])
            if _obj < self._z_int:
                self._z_int = _obj
        else:
            _allint = False
            

        
        return RelaxSolution(obj = _obj,
                             rtime = self._model.RunTime,
                             allint = _allint,
                             pi = _pi,
                             mu = _mu,
                             gamma = 0,
                             omikron_lb = 0,
                             omikron_ub = 0,
                             eta_lb = 0,
                             eta_ub = 0
        )

    
    def log(self,msg):
        if self._model is None:
            return
        self._model.message(msg)

       
    def print_solution_relax(self):
        a = [] 
        for v in self._lbd.select():
            if v.X > .5:
                a.append(v.VarName)

        print(" ".join( a))

            
 
    def column_prepare(self,machine, pi, mu, sigma, gamma):
        nproc = self._instance.nproc
        P = self._instance.P
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S

        x0 = self._instance.mach_map_assign(machine)

        #self._mach[machine].model.reset()

        self._mach[machine]._cte.Obj = - mu
        
        for p in P:
            self._mach[machine]._x[p].Obj = - pi[p]
            if self._args.validate:
                self._mach[machine].model.chgCoeff(
                    self._mach[machine]._pixp_constr,
                      self._mach[machine]._x[p],
                      -pi[p]
                      )


        self._mach[machine].model.update()

        for p in P:
            self._mach[machine]._x[p].Start = x0[p]

    def column_compute(self,machine):
            
        nproc = self._instance.nproc
        P = self._instance.P
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S

        x0 = self._instance.mach_map_assign(machine)

        self._mach[machine].model.optimize(self._cb)
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
            [1* (self._mach[machine]._x[p].X > .5) for p in P], dtype=np.int32
        )
        z = np.array(
            [1* (self._mach[machine]._z[p].X > .5) for p in P], dtype=np.int32
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
            _pixp = self._mach[machine]._pixp.X
            _pmc =self._mach[machine]._pmc.X
            _mmc = self._mach[machine]._mmc.X
            
        return CGColumn(rc=self._mach[machine].model.objVal,
                        rtime = self._mach[machine].model.Runtime,
                        obj=round(self._mach[machine]._obj.X),
                        procs=q,
                        servs=None,
                        g= None,
                        z= None,
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

    def lp_add_col(self,machine,compcol):

        P = self._instance.P
        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        S = self._instance.S
        iL = self._instance.iL
        iN = self._instance.iN

        for v in self._lbd[machine]:
            if np.array_equal(compcol.procs,v._procs):
                return CGAdd(status=CGAddStatus.Exist,var=v)

        
        _col = Column()
        _col.addTerms(
            compcol.procs,
            [self._p_constr[p] for p in P]
            )
        _col.addTerms(1,self._m_constr[machine])

        var = self._model.addVar(
            obj = compcol.obj,
            vtype=GRB.CONTINUOUS,
            ub=1,lb=0,
            column = _col,
            name = "lbd[%d][%d]" % (machine,len(self._lbd[machine]))
        )
        var._procs = compcol.procs
        self._lbd[machine].append(var)
        
        return CGAdd(status=CGAddStatus.Added,var=var)
