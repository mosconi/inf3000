
class Instance:
    def __new__(cls, model=None, assign=None):
        if not model:
            raise ValueError("Missing model")
        if not assign:
            raise ValueError("Missing initial assigment")
        return super(Instance,cls).__new__(cls)
    
    def __init__ (self,model=None,assign=None):


        self._loadmodel(model)
        self._loadassign(assign)

    def _loadmodel(self,_modelfile):

        self.M=[]

        f = open(_modelfile, "r")

        lines = [list(map(int,l.rstrip('\n').split())) for l in f]

        f.close()

        import numpy as np

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


    def _loadassign(self, _assignfile):
        f = open(_assignfile, "r")
        line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
        assign=line[:]


    def validate(self,solution=None):

        if not solution:
            raise("")

 
        import numpy as np
        moved_proc = [n for n,k in enumerate(zip(self.assign,solution)) if k[0]!=k[1]]

        print(moved_proc)

        q=[[] for m in range(self.nmach)]

        q[m].append(np.array([solution[p]==m for p in range(self.nproc)],dtype=np.int32))
        _util = (R*q[m][0].reshape(self.nproc,1)).sum(axis=0) 
#        _util_tr =  
        _obj1= _util - C_bar[m]
        _obj1[_obj1<0]=0
        
        _avail = C[m] - _util
        _obj2=np.array([[bT[r1,r2]*_avail[r1] - _avail[r2] for r2 in range(nres)] for r1 in range(self.nres)],dtype=np.int32)
        _obj2[_obj2<0] = 0
    
        _obj3 = WPMC*RHO*q
        _obj3[[p for p in range(self.nproc) if p in mach_assign[m]]]=0

        _obj5 = MU[solution,m]*q

        _obj = (Wlc*_obj1).sum() + (Wbal*_obj2).sum() + (_obj3).sum() + (WMMC*_obj5).sum()

    def objective(self,assign=None):
        
        if not assign:
            raise("")

        moved_proc = [n for n,k in enumerate(zip(assign,solution)) if k[0]!=k[1]]

        

