import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg3 import CG3
from gurobipy import *

def _cb4(model,where):
    if where == GRB.Callback.MIPSOL:
        nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
        if nodecnt > model._nproc:
            model.terminate()

class CG5(CG3):
        
    def __init__(self,instance,args):
        super().__init__(instance,args)
        self._cb = _cb4

    def lp2mip(self):
        self.__lp2mip()

    def build_lpmodel(self):
        super().build_lpmodel()
        self._lp.Params.Method = 2

    def build_column_model(self,machine):
        super().build_column_model(machine)
        self._mach[machine].model.Params.Method = 2
        self._mach[machine].model._nproc = self._instance.nproc

    def srcs(self):
        self.__srcs()

    def __srcs(self):

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        S = self._instance.S
        N = self._instance.N
        iL = self._instance.iL
        iN = self._instance.iN
        WSMC = self._instance.WSMC

        self._lp.update()

        self._3srcs_constr = tupledict()
        import itertools
        _S = S
        __S = len(_S)
        _total = (nproc)*(nproc-1)*(nproc-2)/(6)
        #_S = {k:v for k,v in S.items() if len(v) > 1}
        _c = 0
        rhs1 = 1
        for i in itertools.combinations(sorted(_S,key=lambda x: len(_S[x]),reverse=True),3):
            for j in itertools.product(S[i[0]],S[i[1]],S[i[2]]):
                _c+=1
                if _c % __S ==0:
                    print("  %15d de %15d" %(_c, _total))

                expr = LinExpr()
                for _lbd in self._lbd.select():
                    expr.addTerms( int((
                        self._lp.getCoeff(self._p_constr[j[0]],_lbd) +
                        self._lp.getCoeff(self._p_constr[j[1]],_lbd) +
                        self._lp.getCoeff(self._p_constr[j[2]],_lbd)
                    )*0.5),
                            _lbd )
                #print(expr <= rhs1 )
                self._3srcs_constr[j] = self._lp.addConstr( expr <= rhs1, name="3src[%d,%d,%d]" % j )

            break
        print("   adicionado %d cortes" % (_c))

        self._lp.update()

    def __srcs2(self):

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        S = self._instance.S
        N = self._instance.N
        iL = self._instance.iL
        iN = self._instance.iN
        WSMC = self._instance.WSMC

        self._lp.update()

        self._3srcs_constr = tupledict()
        import itertools
        _S = S
        __S = len(_S)
        _total = (__S)*(__S-1)*(__S-2)/(6)
        #_S = {k:v for k,v in S.items() if len(v) > 1}
        _c = 0
        for i in itertools.combinations(sorted(_S,key=lambda x: len(_S[x]),reverse=True),3):
            _c+=1
            if _c % __S ==0:
                print("  %15d de %15d" %(_c, _total))
            #print([S[i[0]],S[i[1]],S[i[2]]])
            #print([len(S[i[0]]),len(S[i[1]]),len(S[i[2]])])
            rhs1 = int(len(S[i[0]])/2 + len(S[i[1]])/2 + len(S[i[2]])/2)
            expr = LinExpr()
            for _lbd in self._lbd.select():
                expr.addTerms( int(
                    sum([self._lp.getCoeff(self._p_constr[p],_lbd) for p in S[i[0]]] )/2 +
                    sum([self._lp.getCoeff(self._p_constr[p],_lbd) for p in S[i[1]]] )/2 +
                    sum([self._lp.getCoeff(self._p_constr[p],_lbd) for p in S[i[2]]] )/2 
                    ),
                    _lbd )
            #print(expr <= rhs1 )
            self._3srcs_constr[i] = self._lp.addConstr( expr >= rhs1, name="3src[%d,%d,%d]" % i )

        self._lp.update()

    def __srcs3(self):

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        S = self._instance.S
        N = self._instance.N
        iL = self._instance.iL
        iN = self._instance.iN
        WSMC = self._instance.WSMC

        self._lp.update()

        self._3srcs_constr = tupledict()
        import itertools
        _total = (nserv)*(nmach)
        _S = {k:v for k,v in S.items() if len(v) > 1}
        _c = 0
        for m in range(nmach):
            for s in _S:
                _c+=1
                if _c % min(nmach,nserv) ==0:
                    print("  %15d de %15d" %(_c, _total))
                    #input('Press Enter')
                    
                expr = LinExpr()
                for _lbd in self._lbd.select(m,"*"):
                    _q = 1+ sum([self._lp.getCoeff(self._p_constr[p],_lbd) for p in S[s]])
                    #print("%s: %f" %(_lbd.VarName, _q))
                    expr.addTerms( int(_q * 0.5), _lbd )
                #print(expr <= 1 )
                self._3srcs_constr[m,s] = self._lp.addConstr( expr >= 1, name="3src[%d,%d]" % (m,s) )

        self._lp.update()

        
    def __lp2mip(self):
        self._mip = self._lp

        for v in self._mip.getVars():
            v.vtype = GRB.BINARY


        for m in range(self._instance.nmach):
            self._lbd[m,0].Start=1
    
    def extend(self):
        self.__extend()

    def __extend(self):

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nres = self._instance.nres
        nserv = self._instance.nserv
        sdep = self._instance.sdep
        delta = self._instance.delta
        L = self._instance.L
        S = self._instance.S
        N = self._instance.N
        iL = self._instance.iL
        iN = self._instance.iN
        WSMC = self._instance.WSMC

        x0 = self._instance.map_assign()
            
        self._h=self._lp.addVars(nserv,len(N),vtype=GRB.CONTINUOUS,
                                  lb=0,ub=1,name="h")

        self._o=self._lp.addVars(nserv,len(L),vtype=GRB.CONTINUOUS,
                                  lb=0,ub=1,name="o")
        
        self._lp.update()

        self._h_lb_constr=self._lp.addConstrs((0 - self._h[s,n] >= 0
                                                for s in sorted(S)
                                                for n in N
                                                for m in N[n]),
                                               name="h_lb_constr")
        self._h_ub_constr=self._lp.addConstrs((0 - self._h[s,n] >=0
                                                for s in sorted(S)
                                                for n in N),
                                               name="h_ub_constr")


        self._o_ub_constr=self._lp.addConstrs((-self._o[s,l] + 0 >=0
                                                for s in sorted(S)
                                                for l in L),
                                               name="o_ub_constr")
        self._o_lb_constr=self._lp.addConstrs((-self._o[s,l] + 0  <=0
                                                for s in sorted(S)
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
                             for s in sorted(S)),
                                                   name="spread")

        
        self._lp.update()
        
        # for each column, complete the new constraints

        
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                procs = np.array(v._procs).flatten()
                z = procs - x0[:,m]
                z[z<0] = 0
                for s in sorted(S):
                     # é binário, por causa da restrição de conflito
                    serv = procs[S[s]].sum()
                    if serv >1:
                        print("s: %d, m: %d, sum: %d" % (s,m,serv))
                        print(v.VarName)
                        print(S[s])
                        for p in S[s]:
                            print(procs[p])
                    gserv = z[S[s]].sum()
                    if gserv >1:
                        print("s: %d, m: %d, gsum: %d" % (s,m,gserv))
                        print(S[s])
                    self._lp.chgCoeff( self._h_lb_constr[s,iN[m],m], v, serv)
                    self._lp.chgCoeff( self._h_ub_constr[s,iN[m]], v, serv)
                    self._lp.chgCoeff( self._o_lb_constr[s,iL[m],m], v, serv)
                    self._lp.chgCoeff( self._o_ub_constr[s,iL[m]], v, serv)

        self._lp.update()

    def presolve(self):
        self.__presolve()
        
    def __presolve(self):

        from collections import deque

        self._mip.Params.NumericFocus=1

        lbds = min([self._lbd.select(m,'*') for m in range(self._instance.nmach)],key=len)
        
        queue=deque()
        unqueue=deque()
        allvars = list()
        for m in range(self._instance.nmach):
            allvars.extend(self._lbd.select(m,'*'))
            
        for v in allvars:
            v._obj = v.obj
            v._star= False
            v.obj=0

        
        
        print("-"*80)
        print("  presolve %d" % (m))
        lbds[0]._star=1
        lbds[0].ub=0

        print("  presolve %d allvars" % (m))
        _allvars  = self._lbd.select(m,'*')
        queue.clear()
        
        v = None
        cond = True;
        while cond:
            print("  presolve %d optimize" % (m))
            self._mip.optimize()
            print("  presolve %d status: %d, len: %d" % (m,self._mip.status, len(queue)))
            if self._mip.status == 2:
                for _v in lbds:
                    if _v.X >.5:
                        v = _v
                for _v in allvars:
                    if _v.X >.5:
                        _v._star=True
            else:
                cond = False
                                
                        
            if v:
                print("  presolve %d len(queue): %d" % (m,len(queue)))
                v.ub = 0
                queue.append(v)
            else:
                cond = False

        for v in queue:
            v.ub=1

        lbds[0].ub=1
       
        counter = 0
        for v in allvars:
            v.obj = v._obj
            if v._star:
                counter +=1
            if not v._star :
                v.ub = 0

        print("counter: %d / %d" % (counter,len(allvars)))
        self._smc = self._mip.addVar(name="smc",vtype=GRB.INTEGER,
                                     lb=0,obj=self._instance.WSMC)
        self._g=self._mip.addVars(self._instance.nserv,vtype=GRB.BINARY,
                                  lb=0,ub=1,name="g")

        self._mip.update()
        S = self._instance.S
                
        self._g_constr = self._mip.addConstrs((-self._g[s] == 0 
                                              for s in sorted(S)),
                                              name="g")

        self._smc_constr = self._mip.addConstrs((self._smc >= self._g[s]
                                                for s in sorted(S)),
                                                name="smc")
        self._mip.update()
        x0 = self._instance.map_assign()
        
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                procs = np.array(v._procs).flatten()
                z = procs - x0[:,m]
                z[z<0] = 0
                for s in sorted(S):
                     # é binário, por causa da restrição de conflito
                    gserv = z[S[s]].sum()
                    if gserv >1:
                        print("s: %d, m: %d, gsum: %d" % (s,m,gserv))
                        print(S[s])
                    self._mip.chgCoeff( self._g_constr[s], v, gserv)

        self._mip.update()

        self._mip.Params.NumericFocus=3

    def printrcs(self):
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                print("%15s %9.6f in (%6.3f,%6.3f) %20.3f"%(v.VarName,v.X,v.lb,v.ub,v.RC))

    def printrcs2(self):
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                if v.X > self._args.tol:
                    print("%15s %9.6f in (%6.3f,%6.3f) %20.3f"%(v.VarName,v.X,v.lb,v.ub,v.RC))

    def filter(self):
        cnt = 0
        for m in range(self._instance.nmach):
            for v in self._lbd.select(m,'*'):
                if v.RC> self._args.epslon:
                    if v.ub > self._args.epslon:
                        cnt+=1
                        v.ub=0
                        print("removing: %s %6.3f (%6.3f,%6.3f) %20.3f"%(v.VarName,v.X,v.lb,v.ub,v.RC))
                elif v.ub < self._args.epslon:
                    cnt+=1
                    v.ub = 1
                    print("readding: %s %6.3f (%6.3f,%6.3f) %20.3f"%(v.VarName,v.X,v.lb,v.ub,v.RC))

        self._lp.update()
        return cnt

