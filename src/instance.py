
class Instance:
    def __new__(cls, model=None, assign=None):
        if not model:
            raise ValueError("Missing model")
        if not assign:
            raise ValueError("Missing initial assigment")
        return super(Instance,cls).__new__(cls)
    
    def __init__ (self,model=None,assign=None):
        modelfile=model
        assignfile=assign

        import numpy as np


        self.M=[]

        f = open(modelfile, "r")

        lines = [list(map(int,l.rstrip('\n').split())) for l in f]

        self.nres = lines.pop(0)[0]

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

    
        self.nserv = lines.pop(0)[0]
        delta=[0 for s in range(nserv)]
        sdep={}

        for s in range(self.nserv):
            l = lines.pop(0)
            delta[s]=l[0]
            if l[1]>0:
                sdep[s]=l[2:]

        self.nproc = lines.pop(0)[0]

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
    
        self.R = np.array(R,dtype=np.int32)
        self.S = [s for s in S if s]    

        self.RHO=np.array(RHO,dtype=np.int32)

