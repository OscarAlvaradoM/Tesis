#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Example 8.2 from Malalasekera Book (An Introduction
to Computational Fluid Dynamics THE FINITE VOLUME METHOD 2nd ED)
----------------------------------

A thin plate is initially at a uniform temperature of 200°C. At a certain time
t = 0 the temperature of the east side of the plate is suddenly reduced to 0°C.
The other surface is insulated. Use the explicit finite volume method in con-
junction with a suitable time step size to calculate the transient temperature
distribution of the slab and compare it with the analytical solution at time
(i) t = 40 s, (ii) t = 80 s and (iii) t = 120 s.

In this code case (iii) is solved.

"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np

#-------------Definición de datos iniciales ----------------
lx = 0.02 # meters
flux = 0.0 # °C 
TB = 0. # °C 
k  = 10. # W/m.K
nvx  =  5 # Número de volumenes
rhoC = 10.**7
dt=2.
phi_old=200*np.ones(nvx)
#-------------------------------------------------------------

# ------------Creamos la malla de acuerdo a la geometría definida por el problema
malla = fvm.Mesh(1, volumes=nvx,lengths=lx )
#malla.setDominio(np.array([0., 0.1,0.19, 0.8, 1.2, 1.3, 1.8, 2.0])/100)
malla.info() # Se imprime información

#-------------Definimos fronteras ----------
malla.tagWestWall('N1',flux)
malla.tagEastWall('D1',TB)
#malla.draw()

#--- creamos el objeto coeficientes --------------------------

for i in range(60):
    coef=fvm.Coefficients(malla)
    coef.setDiffusion(k)
    coef.setTime(phi_old,dt,rhoC)
    sistema= fvm.EqSystem(coef)
    sistema.NaiveSetMatrix()
    sol = np.linalg.solve( sistema.mat(),sistema.b() )
    phi_old=sol

print("La sol es:")
print(sol)

graficador = fvm.SolPlotr(malla)
graficador.energy(sol)


