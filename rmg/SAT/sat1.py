import numpy as np

from satispy import Variable,Cnf
from satispy.io import DimacsCnf
from satispy.solver import Minisat

from .sat import SAT

from .satsol import SATSolution

class SAT1(SAT):

    def __init__(self, instance):

        self._instance = instance
        self.global_expr = Cnf()

        #self.solver = Minisat('minisat %s %s')
        self.solver = Minisat()

    def build_model(self):


        nproc = self._instance.nproc
        nmach = self._instance.nmach
        nserv = self._instance.nserv
        nneigh = len(self._instance.N)
        nloc = len(self._instance.L)

        P = self._instance.P
        M = self._instance.M
        S = self._instance.S
        N = self._instance.N
        L = self._instance.L

        io = DimacsCnf()
        
        self.x = np.empty((nproc,nmach), dtype=object)
        for p in P:
            for m in M:
                self.x[p,m]=Variable('x[%d,%d]' % (p,m))
            
        # proc alloc
        all_proc_expr = Cnf()


        for p in P:
            p_expr1 =Cnf()
            p_expr2 =Cnf()
            p_expr3 =Cnf()

            for m1 in M:
                if any(self._instance.R[p] > self._instance.C[m1]):
                    p_expr1 &= -self.x[p,m1]
                else:
                    p_expr1 |= self.x[p,m1]

                    for m2 in range(m1+1,nmach):
                        p_expr2 &= self.x[p,m1] >> - self.x[p,m2]


            all_proc_expr &= p_expr1 & p_expr2 & p_expr3
            

        self.global_expr &= all_proc_expr
        
        # conflict
        all_conflict_expr = Cnf()

        for m in M:
            for s in S:
                if len(S[s]) ==1: continue
                conflict_expr = Cnf()
                for p1 in range(len(S[s])):
                    for p2 in range(p1+1,len(S[s])):
                        conflict_expr &= self.x[S[s][p1],m] >> - self.x[S[s][p2],m]
                all_conflict_expr &= conflict_expr

        self.global_expr &= all_conflict_expr
        
    def solve(self):

        sat_sol = self.solver.solve(self.global_expr)
        
        if not sat_sol.success:
            return SATSolution(False,None)
            
        assignment = np.zeros((self._instance.nproc,self._instance.nmach),dtype=np.int32)
        
        for p in self._instance.P:
            for m in self._instance.M:

                assignment[p,m] = sat_sol[self.x[p,m]]
                
                
        return SATSolution(True,assignment)
            
         
    def exclude_column(self,mach,col):
        expr = Cnf()

        for p in self._instance.P:
            if col[p] ==1 :
                expr |= -self.x[p,mach]
            else:
                expr |=  self.x[p,mach]

        self.global_expr &= expr
                

    def exclude_solution(self,solution):

        expr = Cnf()

        for p in self._instance.P:
            for m in self._instance.P:
                if solution.assignment[p,m] == 1 :
                    expr |= -self.x[p,m]
                else:
                    expr |=  self.x[p,m]
                    
        self.global_expr &= expr
