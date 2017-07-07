from . import *
from gurobipy import *


import numpy as np

class MIP(object):
    _instance = None
    _args = None
    _mip = None
    _lp = None
    _env = None
    def __init__(self,instance,args):
        self._instance = instance
        self._args = args

    def build_model(self):
        raise NotImplementedError("Classe de interface")

    def relax(self):
        if not self._mip:
            raise Exception("MIP model not defined")

        self._mip.relax()
        
    def write(self):
        if not self._mip:
            raise Exception("MIP model not defined")

        self._mip.write("%s.lp" %self._args.run_name)
        self._mip.write("%s.mps" %self._args.run_name)

    def write(self):
        if not self._mip:
            raise Exception("MIP model not defined")

        self._mip.write("%s.lp" %self._args.run_name)
        self._mip.write("%s.mps" %self._args.run_name)
        
    
    
    def build_lp_model(self):
        raise NotImplementedError("Classe de interface")

    def solve(self):
        self._mip.optimize()
        
        _obj = self._mip.ObjVal
        _x = np.array([[self._x[p,m].X > .5 for m in range(self._instance.nmach)] for p in range(self._instance.nproc)],dtype=np.int32)
        _assign = np.empty(self._instance.nproc,dtype=np.int32)
        for p in range(self._instance.nproc):
            for m in range(self._instance.nmach):
                if _x[p,m] == 1:
                    _assign[p] = m
                    continue

        return Solution(obj=_obj, X = _x, assign = _assign)
        

    def solve_lp(self):
        raise NotImplementedError("Classe de interface")
        if not self._lp:
            raise Exception("LP model not defined")
        self._lp.optimize()
        
        _obj = self._lp.ObjVal
        _x = np.array([[self._x[p,m].X for m in range(self._instance.nmach)] for p in range(self._instance.nproc)],dtype=np.float)
        _assign = np.array([ -1 ] * self._instance.nproc,type=np.int32)

        return Solution(obj=_obj, X = _x, assign = _assign)
