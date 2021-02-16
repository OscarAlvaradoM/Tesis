#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Example 5.1 from Malalasekera Book (An Introduction
to Computational Fluid Dynamics THE FINITE VOLUME METHOD 2nd ED)
----------------------------------
A property φ is transported by means of convection and diffusion through
the one-dimensional domain.The boundary conditions are φ 0 = 1 at x = 0 
and φ L = 0 at x = L. Using
five equally spaced cells and the central differencing scheme for convection
and diffusion, calculate the distribution of φ as a function of x for (i) Case 1:
u = 0.1 m/s, (ii) Case 2: u = 2.5 m/s, Case 3: recalculate the solution for 
u = 2.5 m/s with 20 grid nodes and
compare the results with the analytical solution. The following data apply:
length L = 1.0 m, ρ = 1.0 kg/m 3 , Γ = 0.1 kg/m.s.

In this code the Case(iii) is solved
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np

#-------------Definición de datos iniciales ----------------
lx = 1.0 # meters
TA = 1.0 # °C 
TB = 0 # °C 
k  = 0.1 # W/m.K
nvx  =  20 # Número de volumenes
u = 2.5
rho = 1.
#-------------------------------------------------------------

# ------------Creamos la malla de acuerdo a la geometría definida por el problema
malla = fvm.Mesh(1, volumes=nvx,lengths=lx )
#malla.setDominio(np.array([0., 0.1,0.19, 0.8, 1.2, 1.3, 1.8, 2.0])/100)
malla.info() # Se imprime información

#-------------Definimos fronteras ----------
malla.tagWestWall('D1',TA)
malla.tagEastWall('D2',TB)
#malla.draw()

#--- creamos el objeto coeficientes --------------------------
coef=fvm.Coefficients(malla)
coef.setDiffusion(k)
vel = np.array([[u*np.ones(nvx+1) ]])
coef.setAdvection( rho, (vel,) ) #En general 

sistema= fvm.EqSystem(coef)
sistema.NaiveSetMatrix()
print("La matriz que representa el sistema de ecuaciones es:")
print(sistema.mat())

print("EL vector 'b' es:")
print(sistema.b())

sol = np.linalg.solve( sistema.mat(),sistema.b() )
print("La sol es:")
print(sol)

def analytic_sol(x,y=None,z=None):
    return 1. + (1.-np.exp(25*x)) / (7.20*10**10)

#-----Plotting results--------------
graficador = fvm.SolPlotr(malla)
graficador.energy(sol,analytic_sol)