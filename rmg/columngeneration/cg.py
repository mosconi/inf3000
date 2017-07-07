
class CG(object):
    _instance = None
    _args = None
    _mip = None
    _lp = None
    _env = None
    
    def __init__(self,instance,args):
        self._args = args
        self._args = instance

    def build_model(self):
        pass

    def build_lpmodel(self):
        pass
    
    def relax(self):
        pass

    def build_column_model(self,machine):
        pass

    def solve(self):
        pass

    def solve_lp(self):
        pass
