#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Example 8.3 from Malalasekera Book (An Introduction
to Computational Fluid Dynamics THE FINITE VOLUME METHOD 2nd ED)
----------------------------------
Consider convection and diffusion in a one-dimensional domain.
Calculate the transient temperature field if the initial temperature 
is zero everywhere and the boundary conditions are φ = 0 at x = 0 and
∂φ / ∂ x = 0 at x = L. The data are L = 1.5 m, u = 2 m/s, ρ = 1.0 kg/m 3 and
Γ = 0.03 kg/m.s. The source distribution defined by Figure 8.8 applies at
times t > 0 with a = −200, b = 100, x_1 = 0.6 m, x_2 = 0.2 m

Note: In the text the QUICK scheme is used but here the 'upwind1' scheme is
implemented.
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt



#-------------Definición de datos iniciales ----------------
lx = 1.5 # meters
flux = 0.0 # °C 
TA = 0.00 # °C 
k  = 0.03 # W/m.K
nvx  =  45 # Número de volumenes
rho = 1.
u=2.
dt=0.01
phi_old=0.*np.ones(nvx)
vel = np.array([[u*np.ones(nvx+1) ]])
#-------------------------------------------------------------



# ------------Creamos la malla de acuerdo a la geometría definida por el problema
malla = fvm.Mesh(1, volumes=nvx,lengths=lx )
#malla.setDominio(np.array([0., 0.1,0.19, 0.8, 1.2, 1.3, 1.8, 2.0])/100)
#malla.info() # Se imprime información

#-------------Definimos fronteras ----------
malla.tagEastWall('N1',flux)
malla.tagWestWall('D1',TA)
#malla.draw()

def f(x):
    y1 = (x<=0.6)*(-200*x+100)
    y2 = -20 + (20/0.2)*(x-0.6)
    y2 = (0.6<x) * (x<=0.8) * y2
    y3 = (x>0.8) * 0.
    return y1+y2+y3

equis = malla.X()
ye = f(equis)
plt.figure()
plt.plot(equis,ye)

#--- creamos el objeto coeficientes --------------------------


for i in range(300):
    
    coef=fvm.Coefficients(malla)
    coef.setSu(ye)
    coef.setDiffusion(k)
    coef.setAdvection( rho, (vel,), scheme='upwind1')
    coef.setTime(phi_old,dt,rho)
    sistema= fvm.EqSystem(coef)
    sistema.NaiveSetMatrix()
    sol = np.linalg.solve( sistema.mat(),sistema.b() )
    phi_old=sol



print("La sol es:")
print(sol)


graficador = fvm.SolPlotr(malla)
graficador.energy(sol)


