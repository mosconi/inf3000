import numpy as np
from .solution import RelaxSolution,CGColumn,CGValidate,CGValidateStatus,CGAdd,CGAddStatus

from .cg3 import CG3
from gurobipy import *

from itertools import combinations

def _cb4(model,where):
    if where == GRB.Callback.MIPSOL:
        nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
        if nodecnt > model._nproc:
            model.terminate()

class CG5(CG3):
        
    def __init__(self,instance,args):
        super().__init__(instance,args)
        self._cb = _cb4
        nproc = self._instance.nproc
        nmach = self._instance.nmach

        self._ppcounts = np.zeros((nproc, nproc), dtype=np.int32)
        self._mpcounts = np.zeros((nmach, nproc), dtype=np.int32)

    def lp2mip(self):
        self.__lp2mip()

    def build_lpmodel(self):
        super().build_lpmodel()
        self._lp.Params.Method = 2
        #self._lp.Params.Presolve =0 
        self._lp._z_int = np.inf

        nproc = self._instance.nproc
        nmach = self._instance.nmach
        
        for m in range(nmach):
            _lbd = self._lbd[m,0]

            for p0 in range(nproc):
                q_p0 = int(round(self._lp.getCoeff(self._p_constr[p0],_lbd)))
                if q_p0 == 0:
                    continue
                
                self._mpcounts[m,p0]+=q_p0
                for p1 in range(p0+1,nproc):
                    q_p1 = int(round(self._lp.getCoeff(self._p_constr[p1],_lbd)))

                    self._ppcounts[p0,p1]+= (q_p0 + q_p1)//2
                               

            

    def build_column_model(self,machine):
        super().build_column_model(machine)
        self._mach[machine].model.Params.Method = 2
        self._mach[machine].model._nproc = self._instance.nproc

    def cuts_prepare(self):
        nproc = self._instance.nproc
        nmach = self._instance.nmach
        self._3srcs_pp_constr = tupledict()
        self._3srcs_mp_constr = tupledict()

        return 
        for p1 in range(nproc):
            print("%5d PP" % p1)
            for p2 in range(p1+1, nproc):
                c = 0 
                for _lbd in self._lbd.select():
                    if abs(1 - self._lp.getCoeff(self._p_constr[p1],_lbd))<self._args.tol and\
                       abs(1 - self._lp.getCoeff(self._p_constr[p2],_lbd))<self._args.tol:
                            c+=1
                self._ppcounts[p1,p2] = c

            print("%5d MP" % p1)
            for m in range(nmach):
                c=0
                for _lbd in self._lbd.select(m,'*'):
                    if abs(1 - self._lp.getCoeff(self._p_constr[p1],_lbd)) < self._args.tol:
                        c+=1
                self._mpcounts[m,p1] = c

    def cuts_add(self):
        _max_pp = self._ppcounts.max()
        _max_mp = self._mpcounts.max()

        if _max_pp > _max_mp:
            self.__cuts_add_pp()
        else:
            self.__cuts_add_mp()

    def __cuts_add_pp(self):
        nproc = self._instance.nproc
        nmach = self._instance.nmach
        _max_ltp = self._ppcounts.max()
        
        p1,p2 = np.unravel_index(self._ppcounts.argmax(), self._ppcounts.shape)
        print("PP: [%d,%d]: %d" % (p1,p2,self._ppcounts[p1,p2]))
        # p1 é sempre menor que p2, pois é uma matriz triangular superior
        _tp = np.zeros(nproc,dtype=np.int32)
        for p in range(nproc):
            if p == p1: continue
            if p == p2: continue
            if p < p1:
                _tp[p] = self._ppcounts[p1,p2] + self._ppcounts[p,p1] + self._ppcounts[p,p2]
            elif p < p2:
                _tp[p] = self._ppcounts[p1,p2] + self._ppcounts[p1,p] + self._ppcounts[p,p2]
            else:
                _tp[p] = self._ppcounts[p1,p2] + self._ppcounts[p1,p] + self._ppcounts[p2,p]

        _tm = np.zeros(nmach,dtype=np.int32)
        for m in range(nmach):
            _tm[m] = self._ppcounts[p1,p2] + self._mpcounts[m,p1] + self._mpcounts[m,p1]

                
        _max_tp = _tp.max()
        _max_tm = _tm.max()

        import math
        expr = LinExpr()
        if _max_tp > _max_tm:
            _p = _tp.argmax()
            p = tuple(sorted([_p,p1,p2]))
            _min_ltp = min([ self._ppcounts[p[0],p[1]],
                             self._ppcounts[p[0],p[2]],
                             self._ppcounts[p[1],p[2]]
            ])

            for _idx in self._lbd:
                _lbd = self._lbd[_idx]
                print((_idx, _lbd.VarName))
                expr.addTerms(
                    math.floor( (
                        self._lp.getCoeff(self._p_constr[p[0]],_lbd) +
                        self._lp.getCoeff(self._p_constr[p[1]],_lbd) +
                        self._lp.getCoeff(self._p_constr[p[2]],_lbd)
                        ) * 0.5
                    )
                    , _lbd
                )
            self._3srcs_pp_constr[p] = self._lp.addConstr( expr <= 1, name="3src_pp[%d,%d,%d]" % p )
            print("adiconado corte PPP (%d,%d,%d)"%p)
            print("perimetro %d, delta %d" % (_tp.max(),_max_ltp - _min_ltp))
            

        else:
            _m = _tm.argmax()
            p = (_m,p1,p2)
            
            _min_ltp = min([ self._mpcounts[m,p1],
                             self._mpcounts[m,p2]
            ])
            
            for _idx in self._lbd:
                _lbd = self._lbd[_idx]
                print((_idx, _lbd.VarName))
                expr.addTerms(
                    math.floor( (
                        self._lp.getCoeff(self._m_constr[_m],_lbd) +
                        self._lp.getCoeff(self._p_constr[p1],_lbd) +
                        self._lp.getCoeff(self._p_constr[p2],_lbd)
                        ) * 0.5
                    )
                    , _lbd
                )
            self._3srcs_mp_constr[p] = self._lp.addConstr( expr <= 1, name="3src_mp[%d,%d,%d]" % p )
            print("adiconado corte MPP (%d,%d,%d)"%p)
            print("perimetro %d, delta %d" % (_tm.max(),_max_ltp - _min_ltp))

        self._ppcounts[p1,p2] = 0

        self._lp.update()
        
            
        
    def __cuts_add_mp(self):
        nproc = self._instance.nproc
        nmach = self._instance.nmach

        m,p1 = np.unravel_index(self._mpcounts.argmax(), self._mpcounts.shape)
        _max_ltp = self._mpcounts.max()
        print("MP: [%d,%d]: %d" % (m,p1,self._mpcounts[m,p1]))
        
        _tp = np.zeros(nproc,dtype=np.int32)
        for p in range(nproc):
            if p == p1: continue
            if p < p1:
                _tp[p] = self._mpcounts[m,p1] + self._mpcounts[m,p] + self._ppcounts[p,p1]
            else:
                _tp[p] = self._mpcounts[m,p1] + self._mpcounts[m,p] + self._ppcounts[p1,p]

        import math
        expr = LinExpr()

        _p = _tp.argmax()
        if p1 < _p : 
            p = (m,p1,_p)
        else:
            p = (m,_p,p1)

        _min_ltp = min([ self._mpcounts[p[0],p[1]],
                         self._mpcounts[p[0],p[2]],
                         self._ppcounts[p[1],p[2]]
        ])

        self._lp.update()
        self._lbd.clean()
        for _idx in self._lbd:
            print(_idx)
            _lbd = self._lbd[_idx]
            
            print((_idx, _lbd.VarName))
            expr.addTerms(
                math.floor( (
                    self._lp.getCoeff(self._m_constr[m],_lbd) +
                    self._lp.getCoeff(self._p_constr[p1],_lbd) +
                    self._lp.getCoeff(self._p_constr[_p],_lbd)
                ) * 0.5
                )
                , _lbd
            )
        self._3srcs_mp_constr[p] = self._lp.addConstr( expr <= 1, name="3src_mp[%d,%d,%d]" % p )
        print("adiconado corte MPP (%d,%d,%d)"%p)
        print("perimetro %d, delta %d" % (_tp.max(), _max_ltp - _min_ltp))

        self._mpcounts[m,p1] = 0

        self._lp.update()


        
    def cuts_filter(self):
        nproc = self._instance.nproc
        nmach = self._instance.nmach

        zlp =  self._lp.ObjVal
        zinc = self._lp._z_int

        gap  = zinc - zlp
        print(" filter gap %6.3f" % (gap))

        c = 0
        for _idx in self._lbd:
            if _idx[1] == 0 : continue
            # mantém a solução inicial, para evitar problemas infeasible após o
            # corte
            m = _idx[0]
            _lbd = self._lbd[_idx]
            if _lbd.RC > gap and abs(1 - _lbd.ub) < self._args.tol :
                c += 1
                _lbd.ub = 0
                print("removendo %s" %(_lbd.VarName))
                for p1 in range(nproc):
                    qp1 = self._lp.getCoeff(self._p_constr[p1], _lbd)
                    if abs(1 - qp1) < self._args.tol:
                        self._mpcounts[m,p1]-=1
                        for p2 in range(p1,nproc):
                            qp2 =  self._lp.getCoeff(self._p_constr[p2], _lbd)
                            if abs(1 - qp2) < self._args.tol:
                                self._ppcounts[p1,p2]-=1
                self._lp.remove(_lbd)
                                        
        print(" filtered %d vars from %d" %( c , self._lp.NumVars))
        self._lp.update()
        return c
        
    
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


    def lp_add_col(self,mach,colres):
        ares = super().lp_add_col(mach,colres)
        
        nproc = self._instance.nproc

        for p0 in range(nproc):
            q_p0 = colres.procs[p0]
            if q_p0 == 0: continue
            self._mpcounts[mach,p0]+=q_p0
            for p1 in range(p0+1,nproc):
                q_p1 = colres.procs[p1]
                #print((q_p0,q_p1))
                self._ppcounts[p0,p1] += (q_p0 + q_p1 )//2
                #print(self._ppcounts[idx])
            
        for p in range(nproc):
            self._mpcounts[mach,p]+=colres.procs[p]

        return ares
