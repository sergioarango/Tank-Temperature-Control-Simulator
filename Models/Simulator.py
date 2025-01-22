import numpy as np
from Models import ThermoProperties 
from scipy.integrate import odeint
from scipy.optimize import fsolve

class Simulator:
    def __init__ (self, D, H, Tini, t_ini, Q, MW, rho, Patm, Tamb, U):
        self.Thermo = ThermoProperties.ThermoProperties()
        self.D = D
        self.H = H
        self.Patm = Patm
        self.Tamb = Tamb
        self.Tini = Tini
        self.t_ini = t_ini
        self.U = U
        self.Q = [Q, 0]
        self.m = ((np.pi/4)*D**2*H)*rho/MW
        self.A = 2*((np.pi/4)*D**2) + H*np.pi*(D/2)
        self.QA = 0 
        self.resistanceState= 0 

    def Model (self, T_set, tolerance, tp):
        
        def EnergyBalance (H, t, Tin, Q, m, Tamb, U, A):     #Tin: Initial value for temperature,  Q:Resistance power [Watts], m: mass inside tank [mol] 
            dH_dt = Q/m - (U*A*(Tin - Tamb))/m
            return dH_dt
        
        def GetTemperature (var, EnergyBalance, t, Tini, Q, m, Tamb, Patm, U, A):
            T= var
            Hini = self.Thermo.Hliq(T = Tini, P = 0, Patm = Patm)
            H_ii = odeint(EnergyBalance, Hini, t, args=(T, Q, m, Tamb, U, A))
            H_ii = float(H_ii[-1])
            H_s = self.Thermo.Hliq(T = T, P = 0, Patm = Patm)
            solve = H_ii - H_s
            return solve
        
        self.T_set = T_set
        self.tolerance = tolerance
        self.tp = tp
        self.Upper_lim= self.T_set+self.tolerance
        self.Lower_lim= self.T_set-self.tolerance

        t = np.linspace(self.t_ini, self.t_ini + self.tp, 2)
    
        if self.Tini > self.T_set + self.tolerance:
            self.QA = self.Q[1]
            self.QA=float(self.QA)
            self.resistanceState = 0
        elif self.Tini < self.T_set - self.tolerance:
            self.QA = self.Q[0]
            self.QA=float(self.QA)
            self.resistanceState = 1
        
        
        self.current = (self.QA / 110) *1000
        T_system = fsolve(GetTemperature, self.Tini, args=(EnergyBalance, t, self.Tini, self.QA, self.m, self.Tamb, self.Patm, self.U, self.A)) 
        self.Tini = float(T_system[0])
        self.t_ini = self.t_ini + self.tp  
    