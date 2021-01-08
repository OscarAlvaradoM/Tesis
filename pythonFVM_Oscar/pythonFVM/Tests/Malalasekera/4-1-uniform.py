#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Example 4.1 from Malalasekera Book
----------------------------------
Consider the problem of source-free heat conduction in an insulated rod
whose ends are maintained at constant temperatures of 100°C and 500°C
respectively.

The problem solved in the text (Malalasekera) calculates the steady state 
temperature distribution in the rod using a thermal conductivity
k equal to 1000 W/m.K and a cross-sectional area A of 10 × 10−3 m2.
In this code it is assumed that the rod is one dimensional, and therefore it don't 
have a cross-sectional area so a different value of thermal conductivity
(k=10) it's used to compensate. 
In this example a default (uniform except on the borders) mesh is used.

Analytical solution:
    
T = 800 * x + 100
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np

#------------- Initial data definition ----------------
lx = 0.5 # meters
TA = 100 # °C 
TB = 500 # °C 
k  = 10 # W/m.K 
nx  = 10 # number of volumes
#-------------------------------------------------------------

# ------------Mesh definition and border conditions ----
malla = fvm.Mesh(1, volumes=nx,lengths=lx )
malla.info() # Prints some mesh information
malla.tagWestWall('D1',TA) 
malla.tagEastWall('D2',TB)
malla.draw()  #Dislplay a graph of the volumes and nodes in the physical domain

#---- Coefficients procedure --------------------------
coef=fvm.Coefficients(malla)
coef.setDiffusion(k)
#----- System of equations  -----------------------
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
    return 800*x+100

graficador = fvm.SolPlotr(malla)
graficador.energy(sol,analytic_sol)


