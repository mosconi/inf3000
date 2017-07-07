
from gurobipy import Env


class CG(object):
    _instance = None
    _args = None
    _mip = None
    _lp = None
    _env = None
    
    def __init__(self,instance,args):
        self._args = args
        self._instance = instance
        if args.logfile:
            self._env=Env(args.logfile)
        else:
            self._env=Env("columngeneration.log")

        self._mach=dict()
        for m in range(instance.nmach):
            self._mach[m]=lambda: None

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


    def solve(self):
        pass

    def solve_lp(self):
        pass
