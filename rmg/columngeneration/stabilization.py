
import numpy as np
from collections import namedtuple

StabDual = namedtuple('StabDual',["pi",
                                  "mu",
                                  "gamma",
                                  "eta_lb",
                                  "eta_ub",
                                  "omikron_lb",
                                  "omikron_ub"])



class Stabilization(object):

    def __init__(self,instance,args):
        self._instance = instance
        self._args = args
        self._iterations = 0
        self._improvements  = 0
        self._alpha = args.alpha0
        self._best_omega = - np.inf
        self._best_pi = np.zeros(instance.nproc,dtype=np.float64)                        # Duais de processos
        self._best_mu = np.zeros(instance.nmach,dtype=np.float64)                        # Duais de máquina
        self._best_gamma = np.zeros(instance.nserv,dtype=np.float64)                     # Duais de serviço migrado
        self._best_eta_lb = np.zeros((instance.nserv,instance.nmach),dtype=np.float64)       # Dual de h[s,n,m]
        self._best_eta_ub = np.zeros((instance.nserv,len(instance.N)),dtype=np.float64)      # Dual de h[s,n]
        self._best_omikron_lb = np.zeros((instance.nserv,instance.nmach),dtype=np.float64)   # Dual de o[s,l,m]
        self._best_omikron_ub = np.zeros((instance.nserv,len(instance.L)),dtype=np.float64)  # Dual de o[s,l]

    def alpha0(self):
        return self._alpha
    
    def alpha(self,lpres):
        pi =         self._best_pi         + (1 - self._alpha) * lpres.pi
        mu =         self._best_mu         + (1 - self._alpha) * lpres.mu
        gamma =      self._best_gamma      + (1 - self._alpha) * lpres.gamma
        eta_lb =     self._best_eta_lb     + (1 - self._alpha) * lpres.eta_lb
        eta_ub =     self._best_eta_ub     + (1 - self._alpha) * lpres.eta_ub
        omikron_lb = self._best_omikron_lb + (1 - self._alpha) * lpres.omikron_lb
        omikron_ub = self._best_omikron_ub + (1 - self._alpha) * lpres.omikron_ub
        return StabDual(pi,mu,gamma,eta_lb,eta_ub,omikron_lb,omikron_ub)

        
    def linearkp1(self,omega,stabdual):
        
        self._iterations +=1
        new_alpha = max(self._args.alpha_min, self._args.alpha0*(min(1,1 - (self._iterations - self._args.alpha_offset)*(self._args.alpha_scale)/self._instance.nproc)))
        
        if omega > self._best_omega:
            self._improvements += 1
            
            old_alpha = self._alpha
            self._best_pi = stabdual.pi * new_alpha
            self._best_mu = stabdual.mu * new_alpha
            self._best_gamma = stabdual.gamma * new_alpha
            self._best_eta_lb = stabdual.eta_lb * new_alpha
            self._best_eta_ub = stabdual.eta_ub * new_alpha
            self._best_omikron_lb = stabdual.omikron_lb * new_alpha
            self._best_omikron_ub = stabdual.omikron_ub * new_alpha
        
            if self._args.verbose >2:
                print("   omega improment %20.3f -> %20.3f "%( self._best_omega, omega))
            self._best_omega = omega
        self._alpha = new_alpha
            
        return self._alpha

    def linearip1(self,omega,stabdual):
        
        self._iterations +=1
        
        if omega > self._best_omega:
            self._improvements += 1
            new_alpha = max(self._args.alpha_min, self._args.alpha0*(min(1,1 - (self._improvements - self._args.alpha_offset)*(self._args.alpha_scale)/self._instance.nproc)))
            old_alpha = self._alpha
            self._best_pi = stabdual.pi * new_alpha
            self._best_mu = stabdual.mu * new_alpha
            self._best_gamma = stabdual.gamma * new_alpha
            self._best_eta_lb = stabdual.eta_lb * new_alpha
            self._best_eta_ub = stabdual.eta_ub * new_alpha
            self._best_omikron_lb = stabdual.omikron_lb * new_alpha
            self._best_omikron_ub = stabdual.omikron_ub * new_alpha
        
            if self._args.verbose >2:
                print("   omega improment %20.3f -> %20.3f "%( self._best_omega, omega))
                print("   alpha improment %.6f -> %.6f "%( old_alpha, new_alpha))
            self._best_omega = omega
            self._alpha = new_alpha
            
        return self._alpha

    def lineari(self,omega,stabdual):
        
        self._iterations +=1
        
        if omega > self._best_omega:
            self._improvements += 1
            new_alpha = max(self._args.alpha_min, self._args.alpha0*(min(1,1 - (self._improvements - self._args.alpha_offset)*(self._args.alpha_scale)/self._args.alpha_steps)))
            old_alpha = self._alpha
            self._best_pi = stabdual.pi * new_alpha
            self._best_mu = stabdual.mu * new_alpha
            self._best_gamma = stabdual.gamma * new_alpha
            self._best_eta_lb = stabdual.eta_lb * new_alpha
            self._best_eta_ub = stabdual.eta_ub * new_alpha
            self._best_omikron_lb = stabdual.omikron_lb * new_alpha
            self._best_omikron_ub = stabdual.omikron_ub * new_alpha
        
            if self._args.verbose >2:
                print("   omega improment %20.3f -> %20.3f "%( self._best_omega, omega))
                print("   alpha improment %.6f -> %.6f "%( old_alpha, new_alpha))
            self._best_omega = omega
            self._alpha = new_alpha
            
        return self._alpha

    def lineark(self,omega,stabdual):
        
        self._iterations +=1
        new_alpha = max(self._args.alpha_min, self._args.alpha0*(min(1,1 - (self._iterations - self._args.alpha_offset)*(self._args.alpha_scale)/self._args.alpha_steps)))
        
        if omega > self._best_omega:
            self._improvements += 1
            
            old_alpha = self._alpha
            self._best_pi = stabdual.pi * new_alpha
            self._best_mu = stabdual.mu * new_alpha
            self._best_gamma = stabdual.gamma * new_alpha
            self._best_eta_lb = stabdual.eta_lb * new_alpha
            self._best_eta_ub = stabdual.eta_ub * new_alpha
            self._best_omikron_lb = stabdual.omikron_lb * new_alpha
            self._best_omikron_ub = stabdual.omikron_ub * new_alpha
        
            if self._args.verbose >2:
                print("   omega improment %20.3f -> %20.3f "%( self._best_omega, omega))
            self._best_omega = omega
        self._alpha = new_alpha
            
        return self._alpha

    
    def compute(self,omega,stabdual):
        return getattr(self,self._args.method)(omega,stabdual)
    def best_omega(self):
        return self._best_omega
    def improvements(self):
        return self._improvements
