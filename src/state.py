import numpy as np
from model import model as _model

class state:
    def __init__(self, filename: str):
        f=open(filename, "r")
        self.start= [int(i) for i in f.readline().split()]
        self.curr = self.start[:]

    def utilization(self,model: _model, machine: int):
        return model.processes['req'][[p  for p,x in enumerate(self.curr) if x==machine]].sum(axis=0)


    def utilization_all(self,model: _model):
        return np.array([model.processes['req'][[p  for p,x in enumerate(self.curr) if x==m]].sum(axis=0) for m in range(model.nummach)],dtype=int32)

        
