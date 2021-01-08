#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: jose

Example 4.2 from Malalasekera Book (An Introduction
to Computational Fluid Dynamics THE FINITE VOLUME METHOD 2nd ED)
----------------------------------
The problem consists of a large plate of thickness
L = 2 cm with constant thermal conductivity k = 0.5 W/m.K and uniform
heat generation q = 1000 kW/m 3 . The faces A and B are at temperatures
of 100°C and 200°C respectively

In this example the user defines manually the positions of the volumes via the
setDominio() method.
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np

#-------------Initial data ----------------
lx = 0.02 # meters
TA = 100 # °C 
TB = 200 # °C 
k  = 0.5 # W/m.K
nvx  =  5 # Número de volumenes
q = 1000000
#-------------------------------------------------------------

# ------------Mesh creation -----------------
malla = fvm.Mesh(1)
malla.setDominio(np.array([0., 0.1,0.19, 0.8, 1.2, 1.3, 1.8, 2.0])/100)
malla.info() # Se imprime información

#-------------Definition of borders ----------
malla.tagWestWall('D1',TA)
malla.tagEastWall('D2',TB)
malla.draw()
#--- Coefficients procedure--------------------------
coef=fvm.Coefficients(malla)
coef.setSu(q)
coef.setDiffusion(k)

#--------System of equations definition and solving --------
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
    return ((TB-TA)/lx + q*(lx-x)/(2*k))*x +TA

#------- plotting the results-----------
graficador = fvm.SolPlotr(malla)
graficador.energy(sol,analytic_sol)