import argparse

from collections import namedtuple,defaultdict

import numpy as np

parser = argparse.ArgumentParser(add_help=False,
                                 description="Arguments for instance data")

Validate = namedtuple('Validate',['status','obj'])


class Instance:
    def __new__  (self, args):
        if args.instance_filename is None:
            raise Exception("Instance file not provided")
        if args.original_solution_filename is None:
            raise Exception("Originial solution file not provided")

        return super().__new__(self)
    
    def __init__ (self, args):
        
        self.__args=args
        self.__memo=dict()
        self.nproc = 0
        self.nserv = 0
        self.nres = 0
        self.nmach = 0

        self.__loadmodel()
        self.__loadassign()
        
    def __loadmodel(self):

        lines = [list(map(int,l.rstrip('\n').split())) for l in self.__args.instance_filename ]

        nres = lines.pop(0)[0]

        self.nres = nres
        self.SR = range(nres)

        self.T=np.zeros(self.nres,dtype=np.int32)
        self.Wlc=np.zeros(self.nres,dtype=np.int32)


        for r in range(self.nres):
            l = lines.pop(0)
            self.T[r]=l[0]
            self.Wlc[r]=l[1]

        self.nmach = lines.pop(0)[0]

        self.C=np.zeros((self.nmach,self.nres),dtype=np.int32)
        self.SC=np.zeros((self.nmach,self.nres),dtype=np.int32)

        self.MMC=np.zeros((self.nmach,self.nmach),dtype=np.int32)

        
        #L=[[] for i in range(self.nmach)]
        self.L = defaultdict(list)
        self.N = defaultdict(list)

        self.iN = np.array([-1] * self.nmach,dtype=np.int32)
        self.iL = np.array([-1] * self.nmach,dtype=np.int32)
        
        for m in range(self.nmach):
            l = lines.pop(0)
            neigh = l.pop(0) # neighborhood 
            loc = l.pop(0)   # location
            self.iN[m] = neigh
            self.iL[m] = loc
            self.N[ neigh ].append(m) 
            self.L[ loc ].append(m) 
            self.C[m]=l[:self.nres]
            del(l[:self.nres])
            self.SC[m]=l[:self.nres]
            del(l[:self.nres])
            self.MMC[m]=l

        self.M = range(self.nmach)

        nserv = lines.pop(0)[0]
        self.nserv=nserv
        delta=[0 for s in range(nserv)]
        sdep=defaultdict(list)
        sS=defaultdict(list)

        for s in range(self.nserv):
            l = lines.pop(0)
            delta[s]=l[0]
            sS[l[1]].append(s)
            if l[1]>0:
                sdep[s]=l[2:]

        nproc = lines.pop(0)[0]

        self.delta = delta
        self.sdep = sdep
        self.sS = sS
        self.nproc = nproc

        self.S=defaultdict(list)
        self.SM = np.zeros((self.nproc,self.nserv),dtype=np.int32)
        self.R=np.zeros((self.nproc,self.nres),dtype=np.int32)
        self.PMC=np.zeros(self.nproc,dtype=np.int32)
        self.iS=np.zeros(self.nproc,dtype=np.int32)

        
        for p in range(nproc):
            l = lines.pop(0)
            s = l.pop(0)
            self.iS[p] = s 
            self.S[s].append(p)
            self.SM[p,s]=1
            self.R[p]=l[:nres]
            del(l[:nres])
            self.PMC[p]=l[0]

        self.P = range(nproc)
            
        self.nbal = lines.pop(0)[0]
        self.bT=np.zeros((self.nres,self.nres),dtype=np.int32)
        self.Wbal=np.zeros((self.nres,self.nres),dtype=np.int32)
        
        for b in range(self.nbal):
            l = lines.pop(0)
            r1=l[0]
            r2=l[1]
            self.bT[ r1,r2]=l[2]
            l = lines.pop(0)
            self.Wbal[r1,r2] = l[0]

        self.WPMC=lines[0][0]
        self.WSMC=lines[0][1]
        self.WMMC=lines[0][2]
    

    def __loadassign(self):

        line =[list(map(int,l.rstrip('\n').split())) for l in self.__args.original_solution_filename][0]
        self.__assign=line[:]
        self.__mach_assign=[[p[0] for p in enumerate(line) if p[1] == m]  for m in range(self.nmach)]
        self.__mach_map_assign=np.zeros((self.nproc,self.nmach),dtype=np.int32)
        for p in range(self.nproc):
            for m in range(self.nmach):
                self.__mach_map_assign[p,m] = line[p]==m

    def assign(self):
        return self.__assign

    def mach_assign(self,m):
        return self.__mach_assign[m]

    def map_assign(self):
        return self.__mach_map_assign

    def mach_map_assign(self,m):
        return self.__mach_map_assign[:,m]


    def mach_validate(self, machine = None, map_assign = None):

        if machine is None:
            raise Exception("Machine is empty")

        if map_assign is None:
            raise Exception("Map_assign is empty")

        for s in self.S:
            if map_assign[self.S[s]].sum()> 1:
                return Validate(False,0)

        _util = (self.R*map_assign.reshape(self.nproc,1)).sum(axis=0)

        if any(_util > self.C[machine]):
            return Validate(False,0)

        _obj1 = _util - self.SC[machine]
        _obj1[_obj1<0]=0

        _avail = self.C[machine] - _util

        _obj2 = np.array([[
            self.bT[r1,r2]*_avail[r1] - _avail[r2] 
            for r2 in range(self.nres)] 
                          for r1 in range(self.nres)],
                         dtype=np.int32
        )
        _obj2[_obj2<0] = 0
        
        moved_procs = map_assign - self.__mach_map_assign[:,machine]

        moved_procs[moved_procs<0]=0
        
        _obj3=self.PMC*moved_procs
        _obj5=np.array([0],dtype=np.int32)


        _obj=(self.Wlc*_obj1).sum()+(self.Wbal*_obj2).sum()+_obj3.sum()+_obj5.sum()

        return Validate(True,_obj)

    def validate(self,solution=None):

        if not solution:
            raise Exception("")

        if tuple(solution) in self._memo:
            return True
        
        # all proc allocated
        if len(solution) != self.nproc:
            raise Exception("all proc allocated")
            return False

        _sol = np.array(solution,dtype=np.int32)
        moved_proc = np.zeros((self.nmach,self.nproc),dtype=np.int32)

        q=np.zeros((self.nmach,self.nproc),dtype=np.int32)
        q_o=np.zeros((self.nmach,self.nproc),dtype=np.int32)

        for m in range(self.nmach):
            for p in range(self.nproc):
                q[m][p] = _sol[p] ==m
                q_o[m][p] = self.__assign[p] ==m
                if self.__assign[p] == m and solution[p] != m:
                    moved_proc[m][p] =1


        _util = q.dot(self.R)
        _util_tr = moved_proc.dot(self.R*self.T)

        if np.any(_util + _util_tr > self.C):
            raise Exception("usage/transient over capacity")

        
        # conflict constraint + spread constraint
        for s in range(len(self.S)):
            _s = np.array([p in self.S[s] for p in range(self.nproc) ],dtype=np.bool)
            if np.unique(_sol[_s]).size != len(_sol[_s]):
                raise Exception("conflict")
            if len(_sol[_s]) < self.delta[s]:
                raise Exception("spread")


        # dependency

        _n = np.zeros(self.nproc,dtype=np.int32)

        for n in range(len(self.N)):
            _n[np.array([_sol[p] in self.N[n] for p in range(self.nproc)])] = n

        for _s, _dep in self.sdep.items():
            for _p1 in self.S[_s]:
                for _s2 in _dep:
                    if _n[_p1] not in _n[self.S[_s2]]:
                        raise Exception("dependency")
        
        # memoize...
        self._memo[tuple(solution)]={'util': _util,
                                     'moved': moved_proc,
                                     'orig_assign': q_o,
                                     'sol': q
        }

        return True
        

    def objective(self,solution=None):
        
        if not solution:
            raise Exception("")

        if not self.validate(solution):
            raise Exception("")

        _util = self._memo[tuple(solution)]['util']
        q = self._memo[tuple(solution)]['sol']
        q_o = self._memo[tuple(solution)]['orig_assign']
        moved = self._memo[tuple(solution)]['moved']
        
        _obj1= _util - self.SC
        _obj1[_obj1<0]=0

        _avail = self.C - _util

        _obj2=np.zeros((self.nmach,self.nres,self.nres),dtype=np.int32)

        for m in range(self.nmach):
            _obj2[m]=np.array([[self.bT[r1,r2]*_avail[m,r1] - _avail[m,r2] for r2 in range(self.nres)] for r1 in range(self.nres)],dtype=np.int32)
        _obj2[_obj2<0] = 0
        
        _obj3 = self.WPMC*self.PMC*moved
        
        _s = np.zeros(len(self.S),dtype=np.int32)
        _m=moved.sum(axis=0)
        for s in range(len(self.S)):
            _s[s] = _m[self.S[s]].sum()

        
        _obj4=max(_s)

        _obj5 = self.MMC[self.__assign,solution]
        
        _obj = (self.Wlc*_obj1).sum() + (self.Wbal*_obj2).sum() + (_obj3).sum() + (self.WSMC*_obj4) + (self.WMMC*_obj5).sum()

        return _obj
        

    
    def dump(self):
        self.N
