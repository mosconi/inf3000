import numpy as np

def _maptoint(m):
    a = 0
    for b in m:
        a = a*2 + b.item()
    return a

class Instance:
    def __new__(cls, model=None, assign=None,inline=False):
        if not model:
            raise ValueError("Missing model")
        if not assign:
            raise ValueError("Missing initial assigment")
        return super(Instance,cls).__new__(cls)
    
    def __init__ (self,model=None,assign=None,inline=False):

        self._memo={}
        self.__mach_memo = None
        
        self.__inline=inline
        self.__loadmodel(model)
        self.__loadassign(assign)


    def __loadmodel(self,_modelfile):

        if self.__inline:
            lines = [list(map(int,l.rstrip('\n').split())) for l in _modelfile.split('\n')]
        else:
            f = open(_modelfile, "r")
            lines = [list(map(int,l.rstrip('\n').split())) for l in f]
            f.close()

        nres = lines.pop(0)[0]

        self.nres = nres

        self.T=np.zeros(self.nres,dtype=np.bool)
        self.Wlc=np.zeros(self.nres,dtype=np.int32)


        for r in range(self.nres):
            l = lines.pop(0)
            self.T[r]=l[0]
            self.Wlc[r]=l[1]

        self.nmach = lines.pop(0)[0]

        self.C=np.zeros((self.nmach,self.nres),dtype=np.int32)
        self.SC=np.zeros((self.nmach,self.nres),dtype=np.int32)

        self.MU=np.zeros((self.nmach,self.nmach),dtype=np.int32)

        
        L=[[] for i in range(self.nmach)]
        N=[[] for i in range(self.nmach)]

        for m in range(self.nmach):
            l = lines.pop(0)
            N[l.pop(0)].append(m) # neighborhood
            L[l.pop(0)].append(m) # location
            self.C[m]=l[:self.nres]
            del(l[:self.nres])
            self.SC[m]=l[:self.nres]
            del(l[:self.nres])
            self.MU[m]=l

        # remove empty neighborhood/locations
        self.L = [l for l in L if l]
        self.N = [n for n in N if n]

    
        nserv = lines.pop(0)[0]
        self.nserv=nserv
        delta=[0 for s in range(nserv)]
        sdep={}

        for s in range(self.nserv):
            l = lines.pop(0)
            delta[s]=l[0]
            if l[1]>0:
                sdep[s]=l[2:]

        nproc = lines.pop(0)[0]

        self.delta = delta
        self.sdep = sdep
        self.nproc = nproc

        self.S=[]
        self.R=np.zeros((self.nproc,self.nres),dtype=np.int32)
        self.RHO=np.zeros(self.nproc,dtype=np.int32)
        
        S=[[] for i in range(self.nserv)]

        for p in range(nproc):
            l = lines.pop(0)
            S[l.pop(0)].append(p)
            self.R[p]=l[:nres]
            del(l[:nres])
            self.RHO[p]=l[0]
            
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
    
        self.S = [s for s in S if s]    

        self.__mach_memo=[{} for m in range(self.nmach)]


    def __loadassign(self, _assignfile):

        if self.__inline:
            line = [int(l) for l in _assignfile.split()]
        else:
            f = open(_assignfile, "r")
            line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
            f.close()
        self.__assign=line[:]
        self.__mach_assign=[[p[0] for p in enumerate(line) if p[1] == m]  for m in range(self.nmach)]
        self.__mach_map_assign=np.array([[p==m for p in line ] for m in range(self.nmach)],dtype=np.int32)

    def assign(self):
        return self.__assign

    def mach_assign(self,m):
        return self.__mach_assign[m]

    def map_assign(self):
        return self.__mach_map_assign

    def mach_map_assign(self,m):
        return self.__mach_map_assign[m]


    def mach_validate(self, machine = None, map_assign = None):

        if machine is None:
            raise Exception("Machine is empty")

        if map_assign is None:
            raise Exception("Map_assign is empty")

        if _maptoint(map_assign) in self.__mach_memo[machine]:
            print("coluna já calculada")
            return False

        _util = (self.R*map_assign.reshape(self.nproc,1)).sum(axis=0)

        if any(_util > self.C[machine]):
            return False

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
        
#        for s in self.S:
#            if map_assign[s].sum() > 1:
#                return False

        moved_procs = map_assign - self.__mach_map_assign[machine]

        moved_procs[moved_procs<0]=0
        
        _obj3=self.RHO*moved_procs
        _obj5=np.array([0],dtype=np.int32)

        
        self.__mach_memo[machine][_maptoint(map_assign)]={
            'obj': (self.Wlc*_obj1).sum()+(self.Wbal*_obj2).sum()+_obj3.sum()+_obj5.sum(),
            'util': _util,
            'moved_proc': moved_procs
        }

        return True

    def mach_objective(self, machine = None, map_assign = None):
        if machine is None:
            raise Exception("")

        if map_assign is None:
            raise Exception("")

        if not self.mach_validate(machine, map_assign):
            raise Exception("")

        _obj = self.__mach_memo[machine][_maptoint(map_assign)]['obj']

        return _obj




        
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

        #print(moved_proc)

        # capacity constraint + transient usage
        _util = q.dot(self.R)
        #print(_util)
        _util_tr = moved_proc.dot(self.R*self.T)
        #print(_util_tr)
        if np.any(_util + _util_tr > self.C):
            raise Exception("usage/transient")

        
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
        
        _obj3 = self.WPMC*self.RHO*moved
        
        _s = np.zeros(len(self.S),dtype=np.int32)
        _m=moved.sum(axis=0)
        for s in range(len(self.S)):
            _s[s] = _m[self.S[s]].sum()

        
        _obj4=max(_s)

        _obj5 = self.MU[self.__assign,solution]
        
        _obj = (self.Wlc*_obj1).sum() + (self.Wbal*_obj2).sum() + (_obj3).sum() + (self.WSMC*_obj4) + (self.WMMC*_obj5).sum()

        return _obj
        

    
