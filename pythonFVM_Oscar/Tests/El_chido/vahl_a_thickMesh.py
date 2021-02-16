#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jose

Reproduction of the resutls seen in the paper 
'Natural convection of air in a square cavity 
a bench mark solution' by G. De Vahl Davis

In this code a 40x40 mesh is defined and the 
Ra=1e3 case is solved.
"""

import sys
sys.path.append('../../FVM')
import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import gmres
import time
import VisCoFlow as vis
import CheckMethods as chk

plt.close("all")
#-------------Definición de datos iniciales ----------------
lx = 1.0
ly = 1.0
Ra = 1e3
Pr = 0.71
nvx = 40 ; nvy = 40
rho =1.
# ------------Creamos la malla de acuerdo a la geometría definida por el problema---
num_iter=481
#--------malla T -------------------------
malla = fvm.Mesh(2, volumes=(nvx,nvy),lengths=(lx,ly) )
malla.tagNorthWall('Nn1',0.0)
malla.tagSouthWall('Ns2',0.0)
malla.tagEastWall('D1',0.)
malla.tagWestWall('D2',1.)
#malla.draw()
#---------mallla U ------------
mallaU = fvm.UMesh(malla)
mallaU.tagEastWall('D4',0.)
mallaU.tagWestWall('D5',0.)
mallaU.tagNorthWall('D6',0.)
mallaU.tagSouthWall('D7',0.)
#mallaU.draw()
#---------malla V--------------
mallaV = fvm.VMesh(malla)
mallaV.tagEastWall('D4',0.)
mallaV.tagWestWall('D5',0.)
mallaV.tagNorthWall('D6',0.)
mallaV.tagSouthWall('D7',0.)
#mallaV.draw()
#----------malla P ------------------------
mallaP = fvm.Mesh(2, volumes=(nvx,nvy),lengths=(lx,ly) )
mallaP.tagEastWall('D1',0.)
mallaP.tagWestWall('D2',0.)
mallaP.tagNorthWall('D3',0.)
mallaP.tagSouthWall('D4',0.)

#----------Definición de P,U,V--------------------------------------
P = 1.*np.ones((1,nvy,nvx))
#-------------------------------------------------------------------------
#A,alpha = 20.,1. ; pi=np.pi
A,alpha = 2.,1. ; pi=np.pi
X=malla.uStagDef()[0]  ; Y=malla.Y() ; Z = malla.Z()
ye,zeta,equis=np.meshgrid( Y,Z,X) 
U=-A*np.cos(alpha*pi*ye/ly)*np.sin(alpha*pi*equis/lx)
X=malla.X() ; Y=malla.vStagDef()[1] ; Z = malla.Z()
ye,zeta,equis=np.meshgrid( Y,Z,X) 
V=A*np.sin(alpha*pi*ye/ly)*np.cos(alpha*pi*equis/lx) 
vel = (U,V) #---------Tener couidado con la definicion de vel
#-------------------------------------------------------------------------
solU = U[:,:,1:-1]
solV = V[:,1:-1,]
#------------------------------
#     ECUACIÓN DE ENERGÍA

coefT=fvm.Coefficients(malla)
coefT.setDiffusion(1.0)
coefT.setAdvection( 1.0, (U,V),scheme='upwind1' )
sistemaT= fvm.EqSystem(coefT)
sistemaT.NaiveSetCSR()
A = sistemaT.mat()
M2 = spla.spilu(A) ; 
M = spla.LinearOperator(A.shape, M2.solve)
solT, exitCodeT = gmres( A,sistemaT.b(),M=M )
solT = np.reshape(solT,(1,nvy,nvx))    
    

graficador = fvm.SolPlotr(malla)
#graficador.energy(P)
#graficador.energyVel( solT, (solU,solV), (mallaU.X(),mallaV.Y()) )
#solT = 0*solT

t1=time.time()
for i in range(num_iter):
    

    #------------------------------
    #        ECUACIÓN DE MOMENTO
    dP_dx= (P[:,:,:-1]-P[:,:,1:]) / mallaU.dxU()
    #print("dP/dx=",dP_dx)
    dP_dy= (P[:,:-1,:]-P[:,1:,:]) / mallaV.dyV()
    #print("dP/dx=",dP_dy)
    #----------------------------------------
    U_prev = U[:,:,1:-1]
    V_prev = V[:,1:-1,]
    alpha_u=0.5
    alpha_v=0.1
    alpha_p=0.7e-3
    #------- se interpolan las velocidades para encontrar las velocidades sobre la malla original---
    U_interpol = 0.5*(U[:,:,1:] + U[:,:,:-1]) # ESTO SOLO FUNCIONA PARA MALLAS UNIFORMES
    #print("U interpolated for U mesh:",U_interpol)
    V_interpol = 0.5*(V[:,:,1:] + V[:,:,:-1]) # ESTO SOLO FUNCIONA PARA MALLAS UNIFORMES
    #print("V interpolated for U mesh:",V_interpol)

    coefU=fvm.Coefficients(mallaU)
    coefU.setDiffusion(Pr)
    coefU.setAdvection( 1.0, (U_interpol,V_interpol),scheme='upwind1' )
    coefU.setSu(dP_dx)
    coefU.under(alpha_u,U_prev)
    sistemaU= fvm.EqSystem(coefU)
    sistemaU.NaiveSetCSR()
    A = sistemaU.mat()
    M2 = spla.spilu(A) ; 
    M = spla.LinearOperator(A.shape, M2.solve)
    solU, exitCodeU = gmres( sistemaU.mat(),sistemaU.b(), M=M )
    solU = np.reshape(solU,(1,nvy,nvx-1))
    #print("La sol de U (momento) es:")
    #print(solU)   
    
    U_interpol =  0.5*(U[:,1:,:] + U[:,:-1,:]) # ESTO SOLO FUNCIONA PARA MALLAS UNIFORMES
    #print("U interpolated for V mesh:",U_interpol)
    V_interpol = 0.5*(V[:,1:,:] + V[:,:-1,:]) # ESTO SOLO FUNCIONA PARA MALLAS UNIFORMES
    #print("V interpolated for V mesh:",U_interpol)
    coefV=fvm.Coefficients(mallaV)
    coefV.setDiffusion(Pr)
    coefV.setAdvection( 1.0, (U_interpol,V_interpol),scheme='upwind1' )
    bousnsq_term = Pr*Ra*0.5*(solT[:,1:,:]+solT[:,:-1,:])
    coefV.setSu(dP_dy+bousnsq_term)
    coefV.under(alpha_v,V_prev)
    sistemaV= fvm.EqSystem(coefV)
    sistemaV.NaiveSetCSR()   
    A = sistemaV.mat()
    M2 = spla.spilu(A) ; 
    M = spla.LinearOperator(A.shape, M2.solve)
    solV, exitCodeV = gmres(sistemaV.mat(),sistemaV.b(),M=M )
    solV = np.reshape(solV,(1,nvy-1,nvx))
    #print("La sol de V (momento) es:")
    #print(solV)     
    #---------------------------------------------
    #             ECUACIÓN DE CORRECCIÓN A PRESIÓN
    #---------definición de d ------------------
    diJ=mallaU.areasX()/coefU.aP()
    diJ *= alpha_u  #for under relax
    #print("d=",d)
    dIj=mallaV.areasY() / coefV.aP()
    dIj*=alpha_v #for under relax
    
    coefP=fvm.Coefficients(mallaP)
    aE=1e-10*np.ones(malla.nvzyx())
    aE[:,:,:-1] = rho*diJ*mallaU.areasX()
    coefP.setaE(-aE)
    aW=1e-10*np.ones(malla.nvzyx())
    aW[:,:,1:] = rho*diJ*mallaU.areasX()
    coefP.setaW(-aW)
    aN=1e-10*np.ones(malla.nvzyx())
    aN[:,:-1,:] = rho*dIj*mallaV.areasY()
    coefP.setaN(-aN)
    aS=1e-10*np.ones(malla.nvzyx())
    aS[:,1:,:] = rho*dIj*mallaV.areasY()
    coefP.setaS(-aS)
    
    coefP.setaP(aE+aW+aN+aS)
    coefP.saveSu(rho*(U[:,:,:-1]-U[:,:,1:])*mallaP.areasX()+rho*(V[:,:-1,:]-V[:,1:,:])*mallaP.areasY())
    
    sistemaP= fvm.EqSystem(coefP)
    sistemaP.NaiveSetCSR()
    A = sistemaP.mat()
    M2 = spla.spilu(A) ; 
    M = spla.LinearOperator(A.shape, M2.solve)
    solP, exitCodeP = gmres( sistemaP.mat(),sistemaP.b() ,M=M )
    solP =  np.reshape(solP,(1,nvy,nvx))
    #print("La sol de P' es:",solP)
    
    P = P + alpha_p*solP
    #print("La P corregida es:",P)
    
    U[:,:,1:-1] = solU + diJ*(P[:,:,:-1]-P[:,:,1:]) 
    #print("La U corregida es:",U )

    #alpha_v = 0.6
    #V_corrected = solV + dIj*(P[:,-1:,:]-P[:,1:,:])
    #V[:,1:-1,:] = alpha_v*V_corrected +  (1.-alpha_v)*V[:,1:-1,:]
    V[:,1:-1,:] = solV + dIj*(P[:,:-1,:]-P[:,1:,:]) 
    #print("La V corregida es:",V )

    coefT=fvm.Coefficients(malla)
    coefT.setDiffusion(1.0)
    coefT.setAdvection( 1.0, (U,V),scheme='upwind1' )
    sistemaT= fvm.EqSystem(coefT)
    sistemaT.NaiveSetCSR()
    A = sistemaT.mat()
    M2 = spla.spilu(A) ; 
    M = spla.LinearOperator(A.shape, M2.solve)
    solT, exitCodeT = gmres( A,sistemaT.b(),M=M )
    solT = np.reshape(solT,(1,nvy,nvx))  
    
    if (exitCodeT+exitCodeP+exitCodeU+exitCodeV) != 0:
        print("\n===============")
        print("Salida no exitosa en iteración {}".format(i))
        print("Códigos de salida:")
        print("U:{} - V:{} - P:{} - T:{} \n".format(exitCodeU,exitCodeV,exitCodeP,exitCodeT))
        break
    
    if i%int((num_iter-1)/20)==0:
        
#        graficador.energyVel( solT, (solU,solV), (mallaU.X(),mallaV.Y()), size=15 )
        continuidad=np.abs(U[:,:,1:]-U[:,:,:-1]+V[:,1:,:]-V[:,:-1,:])
        avance= ((i+1)/num_iter)*100
        print("{:.1f}% La continuidad es: {}".format( avance ,np.amax(continuidad))  )

t2=time.time()
print("El tiempo fue de: ",t2-t1)

#--------Hacemos velocidades para número de Nusselt-------
U_extend=np.zeros( (3,nvy+2,nvx+2) ) ; V_extend=np.zeros( (3,nvy+2,nvx+2) )
U_extend[1:-1,1:-1,1:-1]=0.5*(U[:,:,1:] + U[:,:,:-1])
V_extend[1:-1,1:-1,1:-1]=0.5*(V[:,1:,:] + V[:,:-1,:])
#-------------------------------------------------------------------------

solT = np.reshape(solT,(1,nvy,nvx)) 
nu,fluxes=chk.nuLine(solT, malla, 'W', (0.0,0.0),ly=1.0,order=2)
print("El número de Nusselt es: ",nu)
nu,fluxes=chk.nuLine(solT, malla, 'Px', (0.5,0.0),ly=1.0,order=2, vel=(U_extend,V_extend))
print("El número de Nusselt es: ",nu)

#nu, nuLocals =chk.nuLine(solT, malla, 'S', (0.0,0.5),lx=1.0,order=2,vel=(U_extend,V_extend))
#print("El número de Nusselt es: ",nu)
#nu, nuLocals =chk.nuLine(solT, malla, 'Py', (0.0,0.5),lx=1.0,order=2,vel=(U_extend,V_extend))
#print("El número de Nusselt es: ",nu)

levels = np.arange(1,10)/10
graficador.Contour( solT, tieneNiv=1, niveles=levels, title="Curvas de nivel Temperatura")


graficador.energyVel( solT, (solU,solV), (mallaU.X(),mallaV.Y()), size=10, colorV="w" , 
                     marcador="o", title="Temperatura y Velocidad", barra=False)
#plt.figure()
vis.colormap2D(malla,solT,'inferno')



