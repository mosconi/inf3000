
import numpy as np
from gurobipy import *

class CG5:
    def __init__(self,
                 instance=None,
                 historysize=10,
                 epslon=0.1
    ):
        if instance is None:
            raise Exception("instance not defined")

        self.__instance = instance
        self._hsz = historysize
        self.__mip = None
        self.__epslon = epslon
        self.__k = 0

    def __build_model(self):
        if self.__mip: 
            return

        self.__mip = Model("mip")
        self.__mip.ModelSense = GRB.MINIMIZE

        # set up MIP vars
        mach_map_assign = self.__instance.map_assign()
        self.__lbd = [[] for m in range(self.__instance.nmach)]
        self.__q = [[] for m in range(self.__instance.nmach)]

        for m in range(self.__instance.nmach):
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
        
    def compute_column(self,machine, pi, alpha):
        return self.__compute_column(machine, pi, alpha)


    def __compute_column(self,machine, pi, alpha):

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
        model.ModelSense = GRB.MINIMIZE
        
        x = model.addVars(nproc, 
                          vtype=GRB.BINARY, 
                          name="x")

        u = model.addVars(nres,
                          lb=0,
                          ub=C,
                          vtype=GRB.INTEGER, 
                          name="u")

        _r = model.addVars(nproc,nres,
                          vtype=GRB.INTEGER, 
                          name="r")


        d = model.addVars(nres,
                          lb=0,
                          ub=C-SC,
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
                _r[p,r] == R[p,r]*x[p] 
                for p in range(nproc)
                for r in range(nres)
            ), 
            name="r"
        )

        model.addConstrs(
            (
                u[r] == quicksum(_r[p,r] for p in range(nproc)) 
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
        model.addConstrs(
            (
                quicksum(x[p] for p in S[s] ) <= 1
                for s in range(len(S)) if len(S[s]) >1
            ),
            name="conflict"
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
            obj3 == WPMC*quicksum(RHO[p]*x[p]
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
        model.setObjective(
            obj1 + obj2 + obj3 + obj5 - pixp - alpha, GRB.MINIMIZE
        )
        model.Params.OutputFlag=0
        model.optimize()

        q = np.array(
            [ 1 * ( x[p].X > .5 ) for p in range(nproc)],
            dtype=np.int32
        )

        return tuple([int(obj1.X + obj2.X + obj3.X + obj5.X), model.objVal, q, model])

    def __relax(self):
        self.__mip.update()

        master = self.__mip.relax()
        p_constr = [master.getConstrByName("p_alloc[%d]" % p) for p in range(self.__instance.nproc)] 
        m_constr = [master.getConstrByName("m_assign[%d]" % m) for m in range(self.__instance.nmach)] 

        return tuple([master, p_constr, m_constr])


    def build_model(self):
        self.__build_model()

    def solve_relax(self,k=None):
        if not self.__mip:
            self.__build_model()

        (master, p_constr, m_constr) = self.__relax()

        if k is not None:
            master.write("relax_%d.lp" % k)
        master.Params.OutputFlag = 0
        master.optimize()
        
        _obj = master.objVal
        
        _pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        _alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        return tuple([_obj, _pi, _alpha])
        
    def __model_pre_optimize(self):
        if not self.__mip:
            self.__build_model()

        (master, p_constr, m_constr) = self.__relax()

        master.Params.OutputFlag = 0
        master.optimize()

        self._last_obj = master.objVal
        
        pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        # print([[master.getCol(v)] for v in master.getVars() if v.X > .5])

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
        master.Params.OutputFlag=0
        master.optimize()
        
        pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)

        # print([[master.getCol(v)] for v in master.getVars() if v.X > .5])

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
            

    def solve_mip(self):
#        self.__mip.write("cg5.lp")
#        self.__mip.Params.OutputFlag=0
        self.__mip.optimize()

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
            
    def mip_add_col(self, machine=None, col=None, obj=None):

        if machine is None:
            raise Exception("machine not defined")

        if col is None:
            raise Exception("column not defined")

        if obj is None:
            raise Exception("objective component not defined")


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

        starting_col = cg.__instance()

    def solve_lr(self,w=[]):
        (lr, p_constr, m_constr) = self.__relax()
        mu=lr.addVars(len(w),obj=w,name="mu")

        lr.Params.OutputFlag = 0
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
        _res = self.__instance.R*q.reshape(nproc,1)
        _util = (_res).sum(axis=0)
        _obj1 = _util - self.__instance.SC[m]
        _obj1[_obj1<0]=0
        
        _avail = self.__instance.C[m] - _util
        _obj2=np.array([[self.__instance.bT[r1,r2]*_avail[r1] - _avail[r2] for r2 in range(nres)] for r1 in range(nres)])
        _obj2[_obj2<0] = 0
    
        _obj3 = self.__instance.WPMC*self.__instance.RHO*q
        _obj3[[p for p in range(nproc) if p in self.__instance.mach_assign(m)]]=0

        _obj5 = self.__instance.WMMC*self.__instance.MU[self.__instance.assign(),m]*q

        _v1 = sum([pi[p]*q[p] for p in range(nproc)])
        _obj = (self.__instance.Wlc*_obj1).sum() + (self.__instance.Wbal*_obj2).sum() + (_obj3).sum() + (_obj5).sum()


        if abs(ObjVal - (_obj - _v1 - alpha)) > epslon:
            print(">>>> cal mismatch!")
            print("  delta:                  %+15.2f" % (ObjVal - (_obj - _v1)-alpha))
            print("  epslon:                 %+15.2f" % epslon )
            print("")

            print("  obj rc calc             %+15.2f"% (_obj - _v1))
            print("  obj rc grb              %+15.2f"% ObjVal)
                        
            print("")
            print("  obj roadef calc         %+15.2f"%_obj)
            print("  obj roadef grb          %+15.2f"% roadef)
            print("  delta:                  %+15.2f" % (roadef - _obj))
            print("")
   
            obj1 = model.getVarByName("obj1")
            print("  obj1 calc               %+15d"%(self.__instance.Wlc*_obj1).sum())
            print("  obj1 grb                %+15.2f"% obj1.X)
            print("  delta:                  %+15.2f" % ((obj1.X - (self.__instance.Wlc*_obj1).sum() )))
            print("")

            if obj1.X - (self.__instance.Wlc*_obj1).sum() > epslon:
                u = [ model.getVarByName("u[0]") for r in range(nres)]
                d = [ model.getVarByName("d[0]") for r in range(nres)]
                print("  CALC _util")
                print(_util)
                print("  GRB  _util")
                print(np.array([ u[r].X  for r in range(nres)], dtype = np.int64))
                print("   deltas:")
                print(np.array([ u[r].X  for r in range(nres)], dtype = np.int64) - _util)

                print("")

                if any(abs(np.array([ u[r].X  for r in range(nres)], dtype = np.int64) - _util)) >0:
                    print("R*x calc")
                    print(_res)
                    print("R*x grb")
                    _r = np.array([[model.getVarByName("r[%d,%d]" %(p,r)).X for r in range(nres)] for p in range(nproc)],dtype=np.int64)
                    print(_r)
                    print("deltas:")
                    print(_res - _r)
                
                print("  SC")
                print(self.__instance.SC[m])
                print("")
                    
                print("  CALC _d")
                print(_obj1)
                print("  GRB  _d")
                print(np.array([ d[r].X  for r in range(nres)], dtype = np.int64))
                print("")
                    
                print("   deltas:")
                print(np.array([ d[r].X  for r in range(nres)], dtype = np.int64) - _obj1)
                print("")


            obj2 = model.getVarByName("obj2")
            print("  obj2 calc               %+15d"%(self.__instance.Wbal*_obj2).sum())
            print("  obj2 grb                %+15.2f"% obj2.X)
            print("  delta:                  %+15.2f" % ((obj2.X - (self.__instance.Wbal*_obj2).sum() )))
            print("")


            print(model.getVars())
            input("Press Enter to continue...")
        
