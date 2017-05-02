
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

    def __build_model(self):
        if self.__mip: 
            return

        self.__mip = Model("mip")
        self.__mip.ModelSense = GRB.MINIMIZE

        # set up MIP vars
        mach_map_assign = self.__instance.mach_map_assign()
        self.__lbd = [[] for m in range(self.__instance.nmach)]

        for m in range(self.__instance.nmach):
            obj = self.__instance.mach_validate(m, mach_map_assign[m])

            self.__lbd[m].append(
                self.__mip.addVar(obj=obj,
                                  vtype=GRB.BINARY,
                                  name="lbd_%d[0]"%m
                )
            )
            
        
        
        self.__mip.update()

        # set up MIP constraints

        ## all process must be assigned
        for p in range(self.__instance.nproc):
            self.__mip.addConstr(
                quicksum(
                    self.__lbd[m][0]*
                    mach_map_assign[m][p] 
                    for m in range(self.__instance.nmach)
                ) == 1,
                name="p_alloc[%d]"%p
            )

        ## all machine must have an allocation
        for m in range(self.__instance.nmach):
            self.__mip.addConstr(
                self.__lbd[m][0] == 1,
                name="m_assign[%d]"%m
            )

        self.__mip.update()
        
        
        
    def __model_pre_optimize(self):
        if not self.__mip:
            self.__build_model()

        master = self.__mip.relax()
        p_constr = [master.getConstrByName("p_alloc[%d]" % p) for p in range(self.__instance.nproc)] 
        m_constr = [master.getConstrByName("m_assign[%d]" % m) for m in range(self.__instance.nmach)] 

        master.Params.OutputFlag = 0 
        master.optimize()

        pi = np.array([c.Pi for c in p_constr], dtype=np.float64)
        alpha = np.array([c.Pi for c in m_constr], dtype=np.float64)
    
        return (master.ObjVal, pi, alpha)


    def __dual_history(self,slack=0.0):
        if not self.__mip:
            self.__build_model()

        _hpi = np.zeros((self._hsz,self.__instance.nproc), dtype=np.float64)
        _halpha = np.zeros((self._hsz,self.__instance.nmach), dtype=np.float64)

        for h in range(self._hsz):
            (obj,pi,alpha) = self.__model_pre_optimize()

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

    def __box_recenter(self):
        for p in range(self.instance.nproc):
            if abs(self._last_pi[p] - self._pi_lb[p]) < self.__epslon:
                self._pi_lb[p]*=(1 - slack)
            if abs(self._last_pi[p] - self._pi_ub[p]) < self.__epslon:
                self._pi_ub[p]*=(1 + slack)

        for m in range(self.instance.nmach):
            if abs(self._last_alpha[m] - self._alpha_lb[m]) < self.__epslon:
                self._alpha_lb[m]*=(1 - slack)
            if abs(self._last_alpha[m] - self._alpha_ub[m]) < self.__epslon:
                self._alpha_ub[m]*=(1 + slack)

    def __rebox(self):
        for _pi in self._last_pi:
            self._pi_lb = self._last_pi * (1 - slack)
            self._pi_ub = self._last_pi * (1 + slack)
            if abs(self._pi_lb)< self.__epslon:
                self._pi_lb = -1
            if abs(self._pi_ub)< self.__epslon:
                self._pi_lb = +1

                
        for _alpha in self._last_alpha:
            self._alpha_lb = self._last_alpha * (1 - slack)
            self._alpha_ub = self._last_alpha * (1 + slack)
            if abs(self._alpha_lb)< self.__epslon:
                self._alpha_lb = -1
            if abs(self._alpha_ub)< self.__epslon:
                self._alpha_lb = +1


            
    def __solve_boxed(self, xi=0.0):
        self.mip.update()
        master = self.mip.relax()

        rlx_proc_constr = [relax_mdl.getConstrByName("p_alloc[%d]" % p ) for p in range(self.instance.nproc)]
        rlx_mach_constr = [relax_mdl.getConstrByName("m_assign[%d]" % m ) for m in range(self.instance.nmach)]

        for p in range(self.instance.nproc):
            col_ub = Column()
            col_lb = Column()

            col_ub.addTerms([+1], [rlx_proc_constr[p]])
            col_lb.addTerms([-1], [rlx_proc_constr[p]])

            w_ub = relax_mdl.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=xi,name="w_ub[%d]"%p,obj=pi_ub[p], column=col_ub)
            w_lb = relax_mdl.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=xi,name="w_lb[%d]"%p,obj=pi_lb[p], column=col_lb)

        for m in range(self.instance.nmach):
            col_ub = Column()
            col_lb = Column()

            col_ub.addTerms([+1], [rlx_mach_constr[m]])
            col_lb.addTerms([-1], [rlx_mach_constr[m]])

        

        pass

    def dual_history(self):
        self.__dual_history()

    def solve(self):
        self.__dual_history()
        pass
        for _xi in self.xis:
            self.__solve_boxed(_xi)
            self.__rebox()
            

    



    

        

    
