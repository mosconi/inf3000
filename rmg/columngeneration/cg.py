
#from collection import defaultdict
from .solution import CGSolution
from gurobipy import Env
import numpy as np

class CG(object):
    _instance = None
    _args = None
    _mip = None
    _lp = None
    _env = None
    _procs = None
    
    def __init__(self,instance,args):
        self._args = args
        self._instance = instance
        self._env=Env(args.logfile)

        self._mach=dict()
        self._procs=dict()
        for m in range(instance.nmach):
            self._mach[m]=lambda: None
            #self._procs[m]= defaultdict(int)

    def __del__(self):
        for m in range(self._instance.nmach):
            del(self._mach[m].model)
            del(self._mach[m])

        if self._lp:
            del(self._lp)
        if self._mip:
            del(self._mip)
        del(self._env)

    def build_model(self):
        pass

    def build_lpmodel(self):
        pass
    
    def relax(self):
        if not self._mip:
            raise Exception("MIP Model not defined")
        if not self._lp:
            self._lp = self._mip.relax()
        pass

    def build_column_model(self,machine):
        pass

    def write(self):
        if not self._mip:
            raise Exception("MIP Model not defined")

        self._mip.write("%s.lp" %self._args.run_name)
        self._mip.write("%s.mps" %self._args.run_name)


    def lpwrite(self):
        if not self._lp:
            raise Exception("LP Model not defined")

        self._lp.write("%s.lp" %self._args.run_name)
        self._lp.write("%s.mps" %self._args.run_name)
        

    def lpiis(self):
        if not self._lp:
            raise Exception("LP Model not defined")
        
        self._lp.computeIIS()
        self._lp.write("%s.ilp" %self._args.run_name)


    def solve(self):
        if not self._mip:
            raise Exception("MIP Model not defined")

        self._mip.optimize()
        _obj = self._mip.ObjVal

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
        
    def solve_lp(self):
        pass
