import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

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
    
    def build_column_model(self,machine):
        pass

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


    def solve(self):
        if not self._model:
            raise Exception("MIP Model not defined")

        self._mmodel.optimize()
        status = self._mip.Status
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

        x = np.zeros((nproc,nmach),dtype=np.int32)
        assign = np.empty(nproc,dtype=np.int32)
        for m in range(nmach):
            v = [ v for v in self._lbd.select(m,'*')  if v.X>0 ][0]
            p = np.array(v._procs).flatten()
            x[:,m] = p
            assign[p==1] = m

        return CGSolution(obj = _obj, X=x, assign = assign)
        
    

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

            
 
