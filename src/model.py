
import numpy as np

class model:
    def __init__(self,filename):
        f = open(filename, "r");

        lines = [ line.rstrip('\n').split() for line in f]
        lines = [list(map(int,l)) for l in lines]

        self.numres = lines.pop(0)[0]
        resources = lines[:self.numres]
        del lines[:self.numres]

        self.nummach = lines.pop(0)[0]
        machines = lines[:self.nummach]
        del lines[:self.nummach]
        
        self.numserv = lines.pop(0)[0]
        services = lines[:self.numserv]
        del lines[:self.numserv]
        
        self.numproc = lines.pop(0)[0]
        processes = lines[:self.numproc]
        del lines[:self.numproc]

        self.numbal = lines.pop(0)[0]
        balances = lines[:self.numbal*2]
        del lines[:self.numbal*2]

        self.wpmc = lines[0][0]
        self.wsmc = lines[0][1]
        self.wmmc = lines[0][2]
        
        _res = np.dtype([('tr',np.bool),('wlc',np.int32)])
        _mach = np.dtype([
            ('neigh',np.int32),
            ('loc',np.int32),
            ('cap',np.int32,(self.numres)),
            ('safecap',np.int32,(self.numres)),
            ('mmc',np.int32,(self.nummach)),
        ])
        #        _serv = np.dtype()
        _proc = np.dtype([
            ('serv', np.int32),
            ('req', np.int32, (self.numres)),
            ('pmc',np.int32)
        ])
        #        _bal
        self.resources = np.array([tuple(l) for l in resources])
        
        self.machines = np.array([tuple([
            l[0],
            l[1],
            tuple(l[2:2+self.numres]),
            tuple(l[2+self.numres:2+2*self.numres]),
            tuple(l[-self.nummach:])
        ]) for l in machines ], dtype=_mach)
        
        self.processes = np.array([tuple([
            l[0],
            tuple(l[1:1+self.numres]),
            l[-1]
        ]) for l in processes ], dtype=_proc )
        print(self.processes)
