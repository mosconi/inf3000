
import numpy as np
from gurobipy import *

_progress_clock = [ '\b|', '\b/' , '\b-', '\b\\']

def _print_backspace():
    print('\b',end='')
    
def _print_prog_clock(cnt):
    print(_progress_clock[cnt%len(_progress_clock)],end='',flush=True)

def _cb(model,where):
    if where == GRB.Callback.SIMPLEX:
        itcnt = model.cbGet(GRB.Callback.SPX_ITRCNT)
        _print_prog_clock(int(itcnt))

class CG5:
    def __init__(self,
                 name=None,
                 instance=None,
                 historysize=10,
                 epslon=0.1
    ):
        if instance is None:
            raise Exception("instance not defined")

        self.__instance = instance
        self._hsz = historysize
        self.__mip = None
        self.__lp = None
        self.__epslon = epslon
        self.__k = 0
        self.__mach_mdl = [None for m in range(instance.nmach)]
        self.__name=name
                
    def __build_lp_model(self):
        if self.__lp:
            return

        self.__lp = Model("lp")
        
        if self.__name:
            self.__lp.Params.LogToConsole = 0
            self.__lp.Params.OutputFlag = 1
            self.__lp.Params.LogFile="lp_%s.log" % self.__name
        else:
            self.__lp.Params.OutputFlag = 0

        self.__lp.Params.MIPGap = self.__epslon
        self.__lp.Params.ScaleFlag = 0
        self.__lp.Params.Quad = 1

        self.__lp.ModelSense = GRB.MINIMIZE

        # set up LP vars
        mach_map_assign = self.__instance.map_assign()
        self.__lp_lbd = [[] for m in range(self.__instance.nmach)]
        self.__lp_q = [[] for m in range(self.__instance.nmach)]
        self.__lp_cols = [set([]) for m in range(self.__instance.nmach)]

        for m in range(self.__instance.nmach):
            self.__lp_cols[m].add(tuple(mach_map_assign[m]))
            obj = self.__instance.mach_objective(m, mach_map_assign[m])
            self.__lp_lbd[m].append(
                self.__lp.addVar(obj=obj,
                                 lb=0,
                                 ub=1,
                                 vtype=GRB.CONTINUOUS,
                                 name="lbd_%d[0]"%m
                )
            )
        self.__lp.update()

        # set up LP constraints

        ## all process must be assigned
        self.__lp_p_alloc={}
        for p in range(self.__instance.nproc):
            self.__lp_p_alloc[p] = self.__lp.addConstr(
                quicksum(
                    self.__lp_lbd[m][0]*
                    mach_map_assign[m][p] 
                    for m in range(self.__instance.nmach)
                ) == 1,
                name="p_alloc[%d]"%p
            )

        ## all machine must have an allocation
        self.__lp_m_assign={}

        for m in range(self.__instance.nmach):
            self.__lp_m_assign[m]=self.__lp.addConstr(
                self.__lp_lbd[m][0] == 1,
                name="m_assign[%d]"%m
            )

        self.__lp.update()

        
    def __build_model(self):
        if self.__mip: 
            return

        self.__mip = Model("mip")
        
        self.__mip.ModelSense = GRB.MINIMIZE

        if self.__name:
            self.__mip.Params.LogToConsole = 0
            self.__mip.Params.OutputFlag = 1
            self.__mip.Params.LogFile="mip_%s.log" % self.__name
        else:
            self.__mip.Params.OutputFlag = 0

        self.__mip.Params.ScaleFlag = 0
        self.__mip.Params.Quad = 1
        
        # set up MIP vars
        mach_map_assign = self.__instance.map_assign()
        self.__lbd = [[] for m in range(self.__instance.nmach)]
        self.__q = [[] for m in range(self.__instance.nmach)]
        self.__cols = [set([]) for m in range(self.__instance.nmach)]

        for m in range(self.__instance.nmach):
            self.__cols[m].add(tuple(mach_map_assign[m]))
            obj = self.__instance.mach_objective(m, mach_map_assign[m])
            self.__lbd[m].append(
                self.__mip.addVar(obj=obj,
                                  vtype=GRB.BINARY,
                                  name="lbd_%d[0]"%m
                )
            )
        
        self.__mip.update()

        # set up MIP constraints

        ## all process must be assigned
        self.__p_alloc={}
        for p in range(self.__instance.nproc):
            self.__p_alloc[p] = self.__mip.addConstr(
                quicksum(
                    self.__lbd[m][0]*
                    mach_map_assign[m][p] 
                    for m in range(self.__instance.nmach)
                ) == 1,
                name="p_alloc[%d]"%p
            )

        ## all machine must have an allocation
        self.__m_assign={}

        for m in range(self.__instance.nmach):
            self.__m_assign[m]=self.__mip.addConstr(
                self.__lbd[m][0] == 1,
                name="m_assign[%d]"%m
            )

        self.__mip.update()

    def __build_column_model(self,machine):
        nres = self.__instance.nres
        nproc = self.__instance.nproc

        pi = np.ones(nproc,dtype=np.int64)
        alpha=1

        R = self.__instance.R
        RHO = self.__instance.RHO
        bT = self.__instance.bT
        T = self.__instance.T
        C = self.__instance.C[machine]
        SC = self.__instance.SC[machine]
        MU = self.__instance.MU[self.__instance.assign(),machine]
        Wlc = self.__instance.Wlc
        Wbal = self.__instance.Wbal
        WPMC = self.__instance.WPMC
        WMMC = self.__instance.WMMC

        model = Model("machine_%d" % machine)

        if self.__name:
            model.Params.LogToConsole = 0
            model.Params.OutputFlag = 1
            model.Params.LogFile="mach_%s_%d.log" % (self.__name,machine)
        else:
            model.Params.OutputFlag = 0

        model.Params.MIPGap = self.__epslon

        model.Params.Quad = 1
        model.Params.ScaleFlag = 0 
            
        model.ModelSense = GRB.MINIMIZE

        roadef = model.addVar(vtype=GRB.INTEGER, 
                          name="roadef")

        x = model.addVars(nproc, 
                          vtype=GRB.BINARY, 
                          name="x")
        
        u = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="u")

        ut = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="ut")

        d = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="d")

        a = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="a")
                          
        b = model.addVars(nres,nres,
                          lb=0,
                          vtype=GRB.INTEGER, 
                          name="b")
        
        obj1 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj1")
        obj2 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj2")
        obj3 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj3")
        obj5 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj5")

        pixp = model.addVar(vtype=GRB.CONTINUOUS, 
                            name="pixp")

        cte = model.addVar(vtype=GRB.BINARY,
                           obj=1,
                           lb=1,
                           ub=1,
                           name="Constant")

        model.update()

        model.addConstrs(
            (
                u[r] == quicksum(R[p,r]*x[p]  for p in range(nproc)) 
                for r in range(nres)
            ), 
            name="util"
        )

        model.addConstrs(
            (
                ut[r] == T[r] * quicksum(
                    R[p,r]*(1-x[p])  
                    for p in range(nproc)
                    if p in self.__instance.mach_assign(machine)
                ) 
                for r in range(nres)
            ), 
            name="util_transient"
        )

        
        model.addConstrs(
            (
                u[r] + ut[r] <= C[r] 
                for r in range(nres)
            ), 
            name="capacity"
        )


        
        model.addConstrs(
            (
                a[r] == C[r] - u[r]  
                for r in range(nres)
            ), 
            name="avail"
        )

        model.addConstrs(
            (
                b[r1,r2] >= bT[r1,r2]*a[r1] - a[r2]
                for r1 in range(nres)
                for r2 in range(nres)
            ),
            name="balance"
        )
                
        model.addConstrs(
            (
                d[r] >= u[r] - SC[r]
                for r in range(nres)
            ),
            name="overload"
        )

        S =  self.__instance.S
        for s in range(len(S)):
            if len(S[s]) >1:
                model.addConstr(
                    (
                        quicksum(x[p] for p in S[s] ) <= 1
                    ),
                    name="conflict[%d]" % s
                )

        _obj1_constr = model.addConstr(
            obj1 == quicksum(Wlc[r]*d[r]
                for r in range(nres)
            ),
            name="obj1"
        )
        _obj2_constr = model.addConstr(
            obj2 == quicksum(Wbal[r1,r2]*b[r1,r2]
                for r1 in range(nres)
                for r2 in range(nres)
            ),
            name="obj2"
        )

        _obj3_constr = model.addConstr(
            obj3 == WPMC*quicksum(RHO[p]*x[p]
            for p in [p for p in range(nproc) 
                      if p not in self.__instance.mach_assign(machine)]
            ),
            name="obj3"
        )

        _obj5_constr = model.addConstr(
            obj5 == WMMC*quicksum(MU[p]*x[p]
            for p in [p for p in range(nproc) 
                      if p not in self.__instance.mach_assign(machine)]
            ),
            name="obj5"
        )

        _pi_constr =  model.addConstr(
            pixp == quicksum(pi[p]*x[p]
                    for p in range(nproc)
            ),
            name="pixp"
        )

        _alpha_constr =  model.addConstr(
            cte == 1, 
            name="Constant"
        )

        model.addConstr(
            roadef == obj1 + obj2 + obj3 + obj5, 
            name="roadef"
        )
        
        model.setObjective(
            roadef - pixp - alpha*cte, GRB.MINIMIZE
        )

        model._x = x

        model._obj1 = obj1
        model._obj2 = obj2
        model._obj3 = obj3        
        model._obj5 = obj5

        model._pi_constr = _pi_constr
        model._alpha_cte = cte

        self.__mach_mdl[machine] = model

        
    def compute_column(self,machine, pi, alpha,k):
        return self.__compute_column2(machine, pi, alpha,k)

    def __compute_column2(self,machine, pi, alpha,k):
        """
        """

        nproc = self.__instance.nproc
        if self.__mach_mdl[machine] is None:
            self.__build_column_model(machine)

        model = self.__mach_mdl[machine]


        for p in range(nproc):
            model.chgCoeff(model._pi_constr,model._x[p],-pi[p])

        model._alpha_cte.setAttr("Obj", -alpha)
        model.update()

        model.reset()
        #model.update()
        model.optimize(_cb)
        _print_backspace()

        q = np.array(
            [ 1 * ( model._x[p].X > .5 ) for p in range(nproc)],
            dtype=np.int32
        )

        roadef = int(model._obj1.X +model._obj2.X +model._obj3.X +model._obj5.X)

        return tuple([roadef, model.objVal, q, model])


    def __compute_column(self,machine, pi, alpha,k):

        nres = self.__instance.nres
        nproc = self.__instance.nproc

        R = self.__instance.R
        RHO = self.__instance.RHO
        bT = self.__instance.bT
        C = self.__instance.C[machine]
        SC = self.__instance.SC[machine]
        MU = self.__instance.MU[self.__instance.assign(),machine]
        Wlc = self.__instance.Wlc
        Wbal = self.__instance.Wbal
        WPMC = self.__instance.WPMC
        WMMC = self.__instance.WMMC

        model = Model("machine_%d" % machine)
        if self.__name:
            model.Params.LogToConsole = 0
            model.Params.OutputFlag = 1
            model.Params.LogFile="mach_%s_%d.log" % (self.__name,machine)
        else:
            model.Params.OutputFlag = 0

        model.Params.MIPGap = self.__epslon        
        model.Params.Quad = 1
        model.Params.ScaleFlag = 0 
        model.Params.NumericFocus = 2
        
        model.ModelSense = GRB.MINIMIZE

        roadef = model.addVar(vtype=GRB.INTEGER, 
                          name="roadef")
        
        x = model.addVars(nproc, 
                          vtype=GRB.BINARY, 
                          name="x")

        u = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="u")

        d = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="d")

        a = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="a")
                          
        b = model.addVars(nres,nres,
                          lb=0,
                          vtype=GRB.INTEGER, 
                          name="b")
        
        obj1 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj1")
        obj2 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj2")
        obj3 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj3")
        obj5 = model.addVar(lb=0,
                            vtype=GRB.INTEGER, 
                            name="obj5")

        pixp = model.addVar(vtype=GRB.CONTINUOUS, 
                            name="pixp")


        model.update()

        model.addConstrs(
            (
                u[r] == quicksum( R[p,r]*x[p] for p in range(nproc)) 
                for r in range(nres)
            ), 
            name="util"
        )
        
        model.addConstrs(
            (
                a[r] == C[r] - u[r]  
                for r in range(nres)
            ), 
            name="avail"
        )

        model.addConstrs(
            (
                b[r1,r2] >= bT[r1,r2]*a[r1] - a[r2]
                for r1 in range(nres)
                for r2 in range(nres)
            ),
            name="balance"
        )
                
        model.addConstrs(
            (
                d[r] >= u[r] - SC[r]
                for r in range(nres)
            ),
            name="overload"
        )

        S =  self.__instance.S
        for s in range(len(S)):
            if len(S[s]) >1:
                model.addConstr(
                    (
                        quicksum(x[p] for p in S[s] ) <= 1
                    ),
                    name="conflict[%d]" % s
                )

        model.addConstr(
            obj1 == quicksum(Wlc[r]*d[r]
                for r in range(nres)
            ),
            name="obj1"
        )
        model.addConstr(
            obj2 == quicksum(Wbal[r1,r2]*b[r1,r2]
                for r1 in range(nres)
                for r2 in range(nres)
            ),
            name="obj2"
        )

        model.addConstr(
            obj3 == WPMC*quicksum(  RHO[p]*x[p]
            for p in [p for p in range(nproc) 
                      if p not in self.__instance.mach_assign(machine)]
            ),
            name="obj3"
        )

        model.addConstr(
            obj5 == WMMC*quicksum(MU[p]*x[p]
            for p in [p for p in range(nproc) 
                      if p not in self.__instance.mach_assign(machine)]
            ),
            name="obj5"
        )

        model.addConstr(
            pixp == quicksum(pi[p]*x[p]
                    for p in range(nproc)
            ),
            name="pixp"
        )

        model.addConstr(
            roadef == obj1 + obj2 + obj3 + obj5, 
            name="roadef"
        )

        model.setObjective(
            roadef - pixp - alpha, GRB.MINIMIZE
        )
        #model.write("machine_%d_v%d.lp" % (machine,k))
        model.optimize()

        q = np.array(
            [ 1 * ( x[p].X > .5 ) for p in range(nproc)],
            dtype=np.int32
        )

        return tuple([int(roadef.X), model.objVal, q, model])

    def __dump_relax(self):
        self.__lp.update()

        master = self.__lp

        
    
    def __relax(self):
        master = None
        if self.__lp:
            self.__lp.update()
            master = self.__lp
        elif self.__mip:
            self.__mip.update()
            master = self.__mip.relax()
        else:
            raise Exception("build_lp_model or build_model must be called")


        p_constr = [master.getConstrByName("p_alloc[%d]" % p) for p in range(self.__instance.nproc)] 
        m_constr = [master.getConstrByName("m_assign[%d]" % m) for m in range(self.__instance.nmach)] 

        return tuple([master, p_constr, m_constr])


    def build_lp_model(self):
        self.__build_lp_model()

    
    def build_model(self):
        self.__build_model()

    def solve_relax(self,k=None):
        
        (master, p_constr, m_constr) = self.__relax()

        if k is not None:
            master.write("relax_%d.lp" % k)

        master.optimize(_cb)
        _print_backspace()

        s = master.Status
        if s == GRB.Status.INF_OR_UNBD or \
           s == GRB.Status.INFEASIBLE or \
           s == GRB.Status.UNBOUNDED or \
           False:
            return tuple([np.nan, [np.nan] , [np.nan] ])
        _obj = master.objVal

        _pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        _alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        return tuple([_obj, _pi, _alpha, master])
        
    def __model_pre_optimize(self):
        if not self.__mip:
            self.__build_model()

        (master, p_constr, m_constr) = self.__relax()

        master.optimize()

        self._last_obj = master.objVal
        
        pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        _skip = [ False ]* self.__instance.nmach

        for m in range(self.__instance.nmach):
            _skip[m] = False
            (obj_roadef, obj,  q) = self.__compute_column(m, pi, alpha[m])
            print("mach %d %+20.3f %+20.3f " % (m, obj_roadef, obj ))
            if obj > -self.__epslon:
                _skip[m] = True
                continue
                
            if not self.__instance.mach_validate(m,q):
                _skip[m] = True
                continue
                
            col = Column()
            col.addTerms(
                q,
                [self.__p_alloc[p] for p in range(self.__instance.nproc)]
            )
            col.addTerms([1],[self.__m_assign[m]])
                         
            self.__lbd[m].append(
                self.__mip.addVar(obj=obj_roadef,
                                  vtype=GRB.BINARY,
                                  column=col,
                                  name="lbd_%d[%d]" % (m, len(self.__lbd[m]))
                )
            )

            self.__mip.update()
        if all(_skip):
            raise Exception("All collumns generated are not valid.")

        return tuple([master.ObjVal, pi, alpha])


    def __dual_history(self,slack=0.0):
        if not self.__mip:
            self.__build_model()

        _hpi = np.zeros((self._hsz,self.__instance.nproc), dtype=np.float64)
        _halpha = np.zeros((self._hsz,self.__instance.nmach), dtype=np.float64)

        for h in range(self._hsz):
            print("dual history %d" %h)
            try:
                (obj,pi,alpha) = self.__model_pre_optimize()
                print(obj)
            except Exception as e:
                print("break")
                print(e)
                break

            
            _hpi[h] = pi
            _halpha[h] = alpha

        self._pi_lb = _hpi.min(axis=0) * (1 - slack)
        self._pi_ub = _hpi.max(axis=0) * (1 + slack)

        self._alpha_lb = _halpha.min(axis=0) * (1 - slack)
        self._alpha_ub = _halpha.max(axis=0) * (1 + slack)

        self._pi_lb[abs(self._pi_lb)<self.__epslon ] = - self._pi_ub[abs(self._pi_lb)<self.__epslon ]
        self._pi_ub[abs(self._pi_ub)<self.__epslon ] = - self._pi_lb[abs(self._pi_ub)<self.__epslon ]

        self._alpha_lb[abs(self._alpha_lb)<self.__epslon ] = - self._alpha_ub[abs(self._alpha_lb)<self.__epslon ]
        self._alpha_ub[abs(self._alpha_ub)<self.__epslon ] = - self._alpha_lb[abs(self._alpha_ub)<self.__epslon ]

        self._pi_lb[abs(self._pi_lb)<self.__epslon ] = - 1
        self._pi_ub[abs(self._pi_ub)<self.__epslon ] = + 1

        self._alpha_lb[abs(self._alpha_lb)<self.__epslon ] = - 1
        self._alpha_ub[abs(self._alpha_ub)<self.__epslon ] = + 1
        

    def __box_recenter(self, slack=0.5):
        ub = self._pi_ub - self._last_pi
        adj = abs(self._pi_ub[abs(ub)< self.__epslon] * (1 - slack))
        if adj.size > 0:
            self._pi_ub[abs(ub)< self.__epslon] += adj
        
        lb = self._pi_lb - self._last_pi
        adj= abs(self._pi_lb[abs(lb)< self.__epslon] * (1 - slack))
        if adj.size > 0:
            self._pi_lb[abs(lb)< self.__epslon] -= adj
        
        ub = self._alpha_ub - self._last_alpha
        adj = abs(self._alpha_ub[abs(ub)< self.__epslon] * (1 - slack))
        if adj.size > 0:
            self._alpha_ub[abs(ub)< self.__epslon] += adj
        
        lb = self._alpha_lb - self._last_alpha
        adj= abs(self._alpha_lb[abs(lb)< self.__epslon] * (1 - slack))
        if adj.size > 0:
            self._alpha_lb[abs(lb)< self.__epslon] -= adj


    def rebox(self,slack=0.1):
        self.__rebox(slack)
        
    def __rebox(self,slack):
        for _pi in self._last_pi:
            self._pi_lb = self._last_pi * (1 - slack)
            self._pi_ub = self._last_pi * (1 + slack)
            if np.any(np.absolute(self._pi_lb)< self.__epslon):
                self._pi_lb[np.absolute(self._pi_lb)< self.__epslon] = -1
            if np.any(np.absolute(self._pi_ub)< self.__epslon):
                self._pi_ub[np.absolute(self._pi_ub)< self.__epslon] = +1

                
        for _alpha in self._last_alpha:
            self._alpha_lb = self._last_alpha * (1 - slack)
            self._alpha_ub = self._last_alpha * (1 + slack)
            if np.any(np.absolute(self._alpha_lb)< self.__epslon):
                self._alpha_lb[np.absolute(self._alpha_lb)< self.__epslon] = -1
            if np.any(np.absolute(self._alpha_ub)< self.__epslon):
                self._alpha_ub[np.absolute(self._alpha_ub)< self.__epslon] = +1


            
    def __solve_boxed(self, xi=0.0):

        (master, p_constr, m_constr) = self.__relax()

        pi_lb = self._pi_lb
        pi_ub = self._pi_ub

        alpha_lb = self._alpha_lb
        alpha_ub = self._alpha_ub

        
        for p in range(self.__instance.nproc):
            col_ub = Column()
            col_lb = Column()

            col_ub.addTerms([+1], [p_constr[p]])
            col_lb.addTerms([-1], [p_constr[p]])

            w_ub = master.addVar(vtype=GRB.CONTINUOUS,
                                 lb=0,
                                 ub=xi,
                                 name="w_ub[%d]"%p,
                                 obj=pi_ub[p],
                                 column=col_ub)
            w_lb = master.addVar(vtype=GRB.CONTINUOUS,
                                 lb=0,
                                 ub=xi,
                                 name="w_lb[%d]"%p,
                                 obj=pi_lb[p],
                                 column=col_lb)

        for m in range(self.__instance.nmach):
            col_ub = Column()
            col_lb = Column()

            col_ub.addTerms([+1], [m_constr[m]])
            col_lb.addTerms([-1], [m_constr[m]])

            v_ub = master.addVar(vtype=GRB.CONTINUOUS,
                                 lb=0,
                                 ub=xi,
                                 name="v_ub[%d]"%p,
                                 obj=alpha_ub[m],
                                 column=col_ub)
            v_lb = master.addVar(vtype=GRB.CONTINUOUS,
                                 lb=0,
                                 ub=xi,
                                 name="v_lb[%d]"%p,
                                 obj=alpha_lb[m],
                                 column=col_lb)
        

        master.update()
        master.optimize()
        
        pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        _skip = [ False ]* self.__instance.nmach

        for m in range(self.__instance.nmach):
            _skip[m] = False
            (obj_roadef, obj,  q) = self.__compute_column(m, pi, alpha[m])
            print("mach %d %+20.3f %+20.3f " % (m, obj_roadef, obj ))
            if obj > -self.__epslon:
                _skip[m] = True
                continue
                
            if not self.__instance.mach_validate(m,q):
                _skip[m] = True
                continue
                
            col = Column()
            col.addTerms(
                q,
                [self.__p_alloc[p] for p in range(self.__instance.nproc)]
            )
            col.addTerms([1],[self.__m_assign[m]])
                         
            self.__lbd[m].append(
                self.__mip.addVar(obj=obj_roadef,
                                  vtype=GRB.BINARY,
                                  column=col,
                                  name="lbd_%d[%d]" % (m, len(self.__lbd[m]))
                )
            )

            self.__mip.update()
        if all(_skip):
            raise Exception("All collumns generated are not valid.")

        return tuple([master.ObjVal, pi, alpha])


    def dual_history(self,hsz=10,slack=0.5):
        self.__hsz= hsz
        self.__dual_history(slack)


    def solve_boxed(self, xi):
        try:
            (obj, pi, alpha) = self.__solve_boxed(xi)
        except Exception as e:
            raise e
        self._last_pi = pi
        self._last_alpha = alpha
        self.__box_recenter()

        

    def solve(self):
        self.__dual_history()
        pass
        for _xi in self.xis:
            self.__solve_boxed(_xi)
            self.__rebox()
            

    def solve_lp2mip(self,file=None):

        self.__mip = self.__lp
        
        for v in  self.__mip.getVars():
            v.vtype = GRB.BINARY

        for m in range(self.__instance.nmach):
            self.__lp_lbd[m][0].Start=1
            
        if file is not None:
            self.__mip.write(file)

        self.__mip.optimize(_cb)
        _print_backspace()

        _obj = self.__mip.ObjVal
        _X = None
        _alloc = None
        for m in range(self.__instance.nmach):
            _lbd = [ _a for _a in self.__lp_lbd[m] if _a.X > .5][0]
            print(_lbd.VarName)
            _col = self.__mip.getCol(_lbd)
            for c_idx in range(len(self.__lp_p_alloc)):
                for c_col in range(_col.size()):
                    if self.__lp_p_alloc[c_idx] == _col.getConstr(c_col):
                        print(c_idx, end=' ', flush=True)
            print()
        
        return tuple([_obj, _X, _alloc])
            

    def solve_mip(self,file=None):
        if file is not None:
            self.__mip.write(file)

        self.__mip.optimize(_cb)
        _print_backspace()

        _obj = self.__mip.ObjVal
        _X = None
        _alloc = None
        for m in range(self.__instance.nmach):
            _lbd = [ _a for _a in self.__lbd[m] if _a.X > .5][0]
            print(_lbd.VarName)
            _col = self.__mip.getCol(_lbd)
            for c_idx in range(len(self.__p_alloc)):
                for c_col in range(_col.size()):
                    if self.__p_alloc[c_idx] == _col.getConstr(c_col):
                        print(c_idx, end=' ', flush=True)
            print()
        
        return tuple([_obj, _X, _alloc])

    def relax(self):
        return self.__relax()

             
    def map_solution(self):
        for m in range(self.__instance.nmach):
            _lbd = [ _a for _a in self.__lbd[m] if _a.X > .5][0]
            print(_lbd.VarName)
            _col = self.__mip.getCol(_lbd)
            for c_idx in range(len(self.__p_alloc)):
                for c_col in range(_col.size()):
                    if self.__p_alloc[c_idx] == _col.getConstr(c_col):
                        print(c_idx, end=' ', flush=True)
            print()

    def lp_add_col(self, machine=None, col=None, obj=None):
        if machine is None:
            raise Exception("machine not defined")

        if col is None:
            raise Exception("column not defined")

        if obj is None:
            raise Exception("objective component not defined")

        if tuple(col) in self.__lp_cols[machine]:
            return 'E'
            for l in self.__lp_lbd[machine]:
                c2 = np.zeros((self.__instance.nproc),dtype=np.int64)
                c = self.__lp.getCol(l)
                for c3 in range(c.size()):
                    for p in range(self.__instance.nproc):
                        if c.getConstr(c3) == self.__lp_p_alloc[p]:
                            c2[p] =1
                if all(c2 == col):
                    print(l.getAttr("VarName"))
                    print("col obj (o/n) %d   %d"  %( l.getAttr("Obj"), obj))
            print("column already defined on machine %d" % machine)
            return


        self.__lp_cols[machine].add(tuple(col))

        _col = Column()
        _col.addTerms(
            col,
            [self.__lp_p_alloc[p] for p in range(self.__instance.nproc)]
        )
        _col.addTerms([1],[self.__lp_m_assign[machine]])
                         
        self.__lp_lbd[machine].append(
            self.__lp.addVar(obj=obj,
                             vtype=GRB.CONTINUOUS,
                             lb=0, ub=1,
                             column=_col,
                             name="lbd_%d[%d]" % (machine, len(self.__lp_lbd[machine]))
            )
        )
        return 'A'
                    
    def mip_add_col(self, machine=None, col=None, obj=None):

        if machine is None:
            raise Exception("machine not defined")

        if col is None:
            raise Exception("column not defined")

        if obj is None:
            raise Exception("objective component not defined")

        if tuple(col) in self.__cols[machine]:
            for l in self.__lbd[machine]:
                c2 = np.zeros((self.__instance.nproc),dtype=np.int64)
                c = self.__mip.getCol(l)
                for c3 in range(c.size()):
                    for p in range(self.__instance.nproc):
                        if c.getConstr(c3) == self.__p_alloc[p]:
                            c2[p] =1
                if all(c2 == col):
                    print(l.getAttr("VarName"))
                    print("col obj (o/n) %d   %d"  %( l.getAttr("Obj"), obj))
            print("column already defined on machine %d" % machine)
            return


        self.__cols[machine].add(tuple(col))

        _col = Column()
        _col.addTerms(
            col,
            [self.__p_alloc[p] for p in range(self.__instance.nproc)]
        )
        _col.addTerms([1],[self.__m_assign[machine]])
                         
        self.__lbd[machine].append(
            self.__mip.addVar(obj=obj,
                              vtype=GRB.BINARY,
                              column=_col,
                              name="lbd_%d[%d]" % (machine, len(self.__lbd[machine]))
            )
        )

    def generate_companion_columns(self,machine=None, col=None,obj=None):

        if machine is None:
            raise Exception("Machine not defined")

        if col is None:
            raise Exception("Column not defined")

        if obj is None:
            obj = self.__instance.mach_objective(machine,col)

        cols = [[] for m in range(self.__instance.nmach)]
        cols[machine].append(np.copy(self.__instance.mach_map_assign(machine)))

        delta_col = col - cols[machine][0]

        p_plus  = delta_col>0
        p_minus = delta_col<0

        print("+ %d" % p_plus.sum())

        from itertools import product
        
        fb = np.array(list(product([0,1], repeat=len(p_minus))),dtype=np.bool)
        print(fb)

        for m in range(self.__instance.nmach):
            if m == machine: continue
            cols[m].append(np.copy(self.__instance.mach_map_assign(m)))
            cols[m][0][p_plus] = 0
            

        print(cols)
        


    def solve_lr(self,w=[]):
        (lr, p_constr, m_constr) = self.__relax()
        mu=lr.addVars(len(w),obj=w,name="mu")

        lr.optimize()
        
        _obj = lr.objVal
        
        _pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        _alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)
        _mu = np.array([mu[m].X for m in range(len(w))], dtype=np.float64)

        return tuple([_obj, _pi, _alpha,_mu])


    def __compute_rc(self,var=None, pi=None, alpha=None):
        if var is None:
            raise Exception("Var not defined")


        rc = var.Obj
        
        col=self.__mip.getCol(var)

        for p in range(self.__instance.nproc):
            rc -= col.getCoef(self.__p_alloc[p])*pi[p]
        
        for m in range(self.__instance.nmach):
            rc -= col.getCoef(self.__m_assign[m])*alpha[m]

        return rc
        

        
    def validate_column(self,q,m,ObjVal,roadef,pi,alpha,epslon,model):
        return self.__validate_column(q,m,ObjVal,roadef,pi,alpha,epslon,model)

    def __validate_column(self,q,m,ObjVal,roadef,pi,alpha,epslon,model):

        nproc = self.__instance.nproc
        nres = self.__instance.nres
        
        R = self.__instance.R
        RHO = self.__instance.RHO
        bT = self.__instance.bT
        C = self.__instance.C[m]
        SC = self.__instance.SC[m]
        MU = self.__instance.MU[self.__instance.assign(),m]
        Wlc = self.__instance.Wlc
        Wbal = self.__instance.Wbal
        WPMC = self.__instance.WPMC
        WMMC = self.__instance.WMMC

        x = q.reshape(nproc,1)

        _res = self.__instance.R*x
        _util = _res.sum(axis=0)
        _d = _util - SC
        _d[_d<0]=0

        _obj1=(Wlc*_d).sum()
        
        _a = C - _util

        _b = np.empty((nres,nres),dtype=np.int64)
        for r1 in range(nres):
            for r2 in range(nres):
                _b[r1,r2] = bT[r1,r2]*_a[r1] - _a[r2]

        _b[_b<0] = 0
        _obj2=(Wbal*_b).sum()

        _rhox =  RHO*(x.reshape(nproc))

        _rhox[[p for p in range(nproc) if p in self.__instance.mach_assign(m)]] =0
        
        _obj3 = WPMC*_rhox.sum()

        _mu = MU*(x.reshape(nproc))
        _mu[[p for p in range(nproc) if p in self.__instance.mach_assign(m)]] =0

        _obj5 = WMMC*_mu.sum()

        _pixp = sum([pi[p]*q[p] for p in range(nproc)])
        _obj = _obj1 + _obj2+ _obj3 + _obj5


        if abs(ObjVal - (_obj - _pixp - alpha)) > epslon:
            print(">>>> cal mismatch!")
            print("  delta:                  %+15.2f" % (ObjVal - (_obj - _pixp- alpha)))
            print("  epslon:                 %+15.2f" % epslon )
            print("")

            print("  obj rc calc             %+15.2f"% (_obj - _pixp -alpha))
            print("  obj rc grb              %+15.2f"% ObjVal)
                        
            print("")
            print("  obj roadef calc      %+15d"% _obj)
            print("  obj roadef grb       %+15d"% roadef)
            print("  delta:                  %+15.2f" % (roadef - _obj))
            print("")

            pixp = model.getVarByName("pixp")

            delta = pixp.X - _pixp
            
            print("  pi*x calc               %+15.2f" % (_pixp) )
            print("  pi*x grb                %+15.2f" % (pixp.X) )
            print("  delta:                  %+15.2f" % delta)
            print("")

            if abs(delta) > epslon:
                print("delta de pi*xp")
   
            obj1 = model.getVarByName("obj1")
            obj2 = model.getVarByName("obj2")
            obj3 = model.getVarByName("obj3")
            obj5 = model.getVarByName("obj5")

            delta = (obj1.X - (self.__instance.Wlc*_obj1).sum() )

            
            print("  obj1 calc            %+15d"%(self.__instance.Wlc*_obj1).sum())
            print("  obj1 grb                %+15.2f"% obj1.X)
            print("  delta:                  %+15.2f" % delta)
            print("")

            if abs(delta) > epslon:
                u = np.array([ model.getVarByName("u[%d]" % r).X for r in range(nres)],dtype=np.int64)
                d = np.array([ model.getVarByName("d[%d]" % r).X for r in range(nres)],dtype=np.int64)
                print("  CALC _util")
                print(_util)
                print("  GRB  _util")
                print(u)
                print("   deltas:")
                print(u - _util)

                print("")

                if any(abs(u - _util)) >0:
                    print("R*x calc")
                    print(_res)
                    print("R*x grb")
                    _r = np.array([[R[p,r] * model.getVarByName("x[%d]" %(p)).X for r in range(nres)] for p in range(nproc)],dtype=np.int64)
                    print(_r)
                    print("deltas:")
                    print(_res - _r)
                
                print("  SC")
                print(self.__instance.SC[m])
                print("")
                    
                print("  CALC _d")
                print(_d)
                print("  GRB  _d")
                print(d)
                print("")
                    
                print("   deltas:")
                print(d - _d)
                print("")


            delta = (obj2.X - (self.__instance.Wbal*_obj2).sum() )
            
            print("  obj2 calc            %+15d"%(self.__instance.Wbal*_obj2).sum())
            print("  obj2 grb             %+15d"% obj2.X)
            print("  delta:                  %+15.2f" % delta)
            print("")
            if abs(delta) > epslon:
                
                u = np.array([ model.getVarByName("u[%d]" % r).X for r in range(nres)], dtype=np.int64)
                a = np.array([ model.getVarByName("a[%d]" % r).X for r in range(nres)], dtype=np.int64)
                b = np.array([ [model.getVarByName("b[%d,%d]" % (r1,r2)).X for r2 in range(nres)] for r1 in range(nres)], dtype=np.int64)

                print("  Wbal")
                print(self.__instance.Wbal)
                print()
                print("  bT")
                print(self.__instance.bT)
                print()
                
                print("  b calc")
                print(_obj2)
                print("  b grb")
                print(b)
                print("  deltas:")
                print(b - _obj2)

                print("  a calc")
                print(_a)
                print("  a grb")
                print(a)
                print("  deltas:")
                deltas = a - _a
                print(deltas)

                if any(abs(deltas)) > epslon:
                    print("  C")
                    print(self.__instance.C[m])
                    print("  u calc")
                    print(_util)
                    print("  u grb")
                    print(u)
                    print("deltas")
                    print(u - _util)

                

            delta = (obj3.X - _obj3.sum() )
            
            print("  obj3 calc            %+15d"% _obj3)
            print("  obj3 grb             %+15d"% obj3.X)
            print("  delta:                  %+15.2f" % delta)
            print("")
            if abs(delta) > epslon:
                _grb_rhox = np.array([ model.getVarByName("_x[%d]" % p).X for p in range(nproc)], dtype=np.int64)
                print("  rho*x calc")
                print(_rhox)
                print("  rho*x grb")
                print(_grb_rhox)
                print("  rho*x delta")
                print(_grb_rhox - _rhox)
                
                
            delta = (obj5.X - _obj5.sum() )
            
            print("  obj5 calc            %+15d"% _obj5.sum())
            print("  obj5 grb                %+15.2f"% obj5.X)
            print("  delta:                  %+15.2f" % delta)
            print("")
            if abs(delta) > epslon:
                print("delta de obj5")
                

#            model.write("modelo_mismatch.lp")
#            model.write("modelo_mismatch.mps")
#            input("Press Enter to continue...")
#           raise Exception("Mismatch")
        
