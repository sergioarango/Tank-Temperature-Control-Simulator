# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:16:10 2025

@author: sarango
"""

from Models import Simulator
import matplotlib.pyplot as plt

Din = 0.306             #metros
H = 0.1631              #metros
U = 0.49                #W/m2K
Tini = 15 + 273.15     #K
t_ini = 0               #seconds
Q = 0.5                 #kW
MW = 18
rho = 1000              #kg/m3
Patm = 75.24            #Kpa
Tamb = 15 + 273.15      #K
U = 0.05                #kW/m2K
timesimulation = 3600*2 #seconds
T_set = 35 + 273.15     #K
tolerance = 5
tp = 1                  #seconds

tank = Simulator.Simulator(D = Din, H = H, 
                           Tini = Tini, t_ini = t_ini, 
                           Q = Q, MW = MW, rho = rho, Patm = Patm, Tamb = Tamb, U = U)
t = []
T_model = []
while tank.t_ini < timesimulation:
    t.append(tank.t_ini)
    T_model.append(tank.Tini)
    tank.Model(T_set = T_set, tolerance = tolerance, tp = tp)

plt.plot(t, T_model)
plt.show()    
