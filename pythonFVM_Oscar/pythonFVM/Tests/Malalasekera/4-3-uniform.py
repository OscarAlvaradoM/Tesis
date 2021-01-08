#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Example 4.3 from Malalasekera Book (An Introduction
to Computational Fluid Dynamics THE FINITE VOLUME METHOD 2nd ED)
----------------------------------
The problem is a cylindrical fin. The base is at a temperature of 100°C (TB) and 
the end is insulated. The fin is exposed to an ambient temperature of 20°C.
For more details refer to the cited text.

In this example a default (uniform except on the borders) mesh is used.
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np

#------------ Initial data  ----------------
lx = 1. # meters
TA = 100 # °C 
dB = 0. 
T_inf = 20
k  = 1 # W/m.K
nvx  =  5 # Número de volumenes
n=5
#-------------------------------------------------------------

# -----------Mesh definition --------------
malla = fvm.Mesh(1, volumes=nvx,lengths=lx )
malla.info() # Se imprime información
#- Borders definition
malla.tagWestWall('D1',TA)
malla.tagEastWall('N1',dB)
malla.draw()

#--- Coefficients procedure --------------------------
coef=fvm.Coefficients(malla)
coef.setSp(-n**2)
coef.setSu(T_inf*n**2)
coef.setDiffusion(k)
#--------- Section for getting the System of Equations and it's solution ------

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
    return np.cosh(n*(lx-x))*(TA-T_inf)/np.cosh(n*lx)+T_inf

#plotting the solution
graficador = fvm.SolPlotr(malla)
graficador.energy(sol,analytic_sol)

