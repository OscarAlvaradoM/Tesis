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
In this example the user defines manually the positions of the volumes via the
setDominio() method.

Analytical solution:
    
T = 800 * x + 100
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt

#-------------Definición de datos iniciales ----------------
lx = 0.5 # meters
TA = 100 # °C 
TB = 500 # °C 
k  = 10 # W/m.K
nx  = 5 # Número de fronteras (5 volumenes tienen 6 fronteras )
#-------------------------------------------------------------

# ------------Creamos la malla de acuerdo a la geometría definida por el problema
malla = fvm.Mesh(1)
malla.setDominio(np.array([0., 0.05, 0.1, 0.3, 0.48, 0.5]))
malla.info() # Se imprime información

#-------------Definimos fronteras ----------
malla.tagWestWall('D1',TA)
malla.tagEastWall('D2',TB)
print(malla.X(),malla.Y(),malla.Z())
malla.draw()

#--- creamos el objeto coeficientes --------------------------
coef=fvm.Coefficients(malla)
coef.setDiffusion(k)
print("\naP es:")
print(coef.aP())
print("\naE es:")
print(coef.aE())
print("\naW es:")
print(coef.aW())

#---- imprimimos diccionarios con los valores de cada frontera
print("Los valores Dirichlet son:")
print(malla.dirichValues() )


#print("Los índices de los nodos que tienen incógnita son:")
#print(malla.unkown_mesh_indx())

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


#malla.printCoordTags()
#print(malla.deltaX() )

#
# Calculamos la solución analítica
#
#Ta = 800 * x + 100
#
#  Se grafica la solución
#
#x *= 100 # Transformación a [cm]

#---------------Código de Graficcación------------------------------
#plt.plot(x,Ta, '-', label = 'Sol. analítica') # Sol. analítica
#plt.plot(x,T,'o', label = 'Sol. FVM')
#plt.title('Solución de $k (\partial^2 T/\partial x^2) = 0$ con FVM')
#plt.xlabel('$x$ [cm]')
#plt.ylabel('Temperatura [$^o$C]')
#plt.grid()
#plt.legend()
#plt.savefig('example01.pdf')
#plt.show()
