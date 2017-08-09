import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg import CG
from gurobipy import *


class CG4(object):
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

    def __del__(self):
        for m in range(self._instance.nmach):
            if hasattr(self._mach[m],'model'):
                del(self._mach[m].model)
            del(self._mach[m])

        if self._lp:
            del(self._lp)
        if self._mip:
            del(self._mip)
        del(self._env)

    def build_lpmodel(self,ndeps):
        return self.__build_lpmodel(ndeps)

    def __build_lpmodel(self,ndeps):
        
        self._lp = Model("RM_%d" % ndeps, env=self._env)
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
        sS = self._instance.sS
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
        
        xstart = None
        __this_procs = list()

        print("# deps ", ndeps)
        for s in sS[ndeps]:
            print("serv ", s, sdep[s])
            __this_procs += S[s]

        print(__this_procs)
