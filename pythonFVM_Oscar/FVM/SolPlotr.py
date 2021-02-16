#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 16:47:33 2019

@author: jose
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import VisCoFlow as vis

class SolPlotr():
    
    def __init__(self, mesh):
        
        self._mesh = mesh
        self._lx = mesh.lx()
        self._ly = mesh.ly()
        self._lz = mesh.lz()

    def scalarByCut(self,sol,x_cut=None, y_cut=None, z_cut=None,title='Solution',size=20, min_val=None,max_val=None):
        
        mesh=self._mesh
        tags = mesh.tags()
        Nx=mesh.NdomX() ; Ny=mesh.NdomY() ; Nz=mesh.NdomZ()
        tagsCut=np.reshape(tags,(Nz,Ny,Nx))[1:-1,1:-1,1:-1]
        volX=mesh.X() ; volY=mesh.Y() ; volZ=mesh.Z()
        coordY,coordZ,coordX = np.meshgrid(volY,volZ,volX)
        if x_cut==None: x_cut=mesh.lx()/2.
        if y_cut==None: y_cut=mesh.ly()/2.
        if z_cut==None: z_cut=mesh.lz()/2.
        indx_xcut=np.argmax(x_cut<volX)
        indx_ycut=np.argmax(y_cut<volY)
        indx_zcut=np.argmax(z_cut<volZ)
        if (x_cut-volX[indx_xcut-1]) < (volX[indx_xcut] - x_cut): 
            indx_xcut-=1  #Si está mas cerca de X[índice anterior] tomamos el índice anterior
        if (y_cut-volY[indx_ycut-1]) < (volY[indx_ycut] - y_cut): 
            indx_ycut-=1
        if (z_cut-volZ[indx_zcut-1]) < (volZ[indx_zcut] - z_cut): 
            indx_zcut-=1
        #========================================================
        # GRAFICAMOS EL CORTE Y
        #========================================================
        X=coordX[:,indx_ycut,:] ; Y=coordY[:,indx_ycut,:] ; Z=coordZ[:,indx_ycut,:] 
        tags2ndCut = tagsCut[:,indx_ycut,:] ; solCut = sol[:,indx_ycut,:] 
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tags2ndCut.flatten(),'phi':solCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]       
        plt.figure()
        plt.scatter(X,Z,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val)
        plt.title( title+" Plano XZ con y={:.3f}".format(volY[indx_ycut]) )
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.lz())
        plt.colorbar()
        plt.show()
        #========================================================
        # GRAFICAMOS EL CORTE X
        #========================================================
        X=coordX[:,:,indx_xcut] ; Y=coordY[:,:,indx_xcut] ; Z=coordZ[:,:,indx_xcut] 
        tags2ndCut = tagsCut[:,:,indx_xcut] ; solCut = sol[:,:,indx_xcut]
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tags2ndCut.flatten(),'phi':solCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]       
        plt.figure()
        plt.scatter(Y,Z,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val)
        plt.title( title+" Plano YZ con x={:.3f}".format(volX[indx_xcut]) )
        plt.xlim(0.,mesh.ly())
        plt.ylim(0.,mesh.lz())
        plt.colorbar()
        plt.show()
        #========================================================
        # GRAFICAMOS EL CORTE Z
        #========================================================
        X=coordX[indx_zcut,:,:] ; Y=coordY[indx_zcut,:,:] ; Z=coordZ[indx_zcut,:,:] 
        tags2ndCut = tagsCut[indx_zcut,:,:] ; solCut = sol[indx_zcut,:,:]
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tags2ndCut.flatten(),'phi':solCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]       
        plt.figure()
        plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val)
        plt.title( title+" Plano XY con z={:.3f}".format(volZ[indx_zcut]) )
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        plt.colorbar()
        plt.show()

        
    def energy(self,sol,f_anlytc=None,title='Solution',size=20):
        
        mesh=self._mesh
        sol = sol.flatten()
        tags = mesh.tags()
        volX=mesh.X()
        volY=mesh.Y()
        volZ=mesh.Z()
        coordY,coordZ,coordX = np.meshgrid(volY,volZ,volX)
        X=coordX.flatten()
        Y=coordY.flatten()
        Z=coordZ.flatten()
        dim = mesh.dim()
        Nx=mesh.NdomX() ; Ny=mesh.NdomY() ; Nz=mesh.NdomZ()
        tagsCut=np.reshape(tags,(Nz,Ny,Nx))[1:-1,1:-1,1:-1]
        tagsCut=tagsCut.flatten()
        data = {'X':X,'Y':Y,'Z':Z,'Tag':tagsCut,'phi':sol}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        ##graficar primero probar qué está imprimiendo
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]
        
        if f_anlytc:
            X=list(X)
            X=np.array(X)
            exac_phi=f_anlytc(X,Y,Z)
            error = phi-exac_phi
            max_e =np.amax(error)
            distance = (np.dot(error,error))**0.5
        #bool_cond=df['Tag'].str.contains('D', regex=False)
        #dirich_coordTags=np.transpose(df[bool_cond].values)
        #bool_cond=df['Tag'].str.contains('N', regex=False)
        #neum_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------
        if dim ==1:
            plt.figure()
            plt.plot(X,phi,'o',markersize=5)
            plt.title(title)
            plt.xlabel('X')
            plt.ylabel('f(X)')
            plt.xlim(0.,mesh.lx())
            if f_anlytc:
                plt.figure()
                plt.plot(X, 100*error/exac_phi ,'s')
                plt.title('Max Err={:.5f} Dist={:.5f}'.format(max_e,distance))
                plt.xlabel('X')
                plt.ylabel('% Error')
        elif dim ==2:
            plt.figure()
            plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size)
            plt.title(title)
            plt.xlim(0.,mesh.lx())
            plt.ylim(0.,mesh.ly())
            plt.colorbar()
        elif dim ==3:
            fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.scatter(dirich_coordTags[0], dirich_coordTags[1], dirich_coordTags[2], color='y')
            #ax.scatter(neum_coordTags[0],neum_coordTags[1] ,neum_coordTags[2] , color='b')
            #ax.scatter(source_coordTags[0],source_coordTags[1] , source_coordTags[2], color='m')
            #ax.scatter(inter_coordTags[0],inter_coordTags[1] , inter_coordTags[2], color='k')
        plt.show()

    def energyContour(self,sol,f_anlytc=None,title='Solution'):
        
        mesh=self._mesh
        tags = mesh.tags()
        volX=mesh.X()
        volY=mesh.Y()
        volZ=mesh.Z()
        nvx =mesh.nvx() ; nvy = mesh.nvy()
        coordY,coordZ,coordX = np.meshgrid(volY,volZ,volX)
        sol=sol.flatten()
        X=coordX.flatten()
        Y=coordY.flatten()
        Z=coordZ.flatten()
        plt.figure()
        plt.contour(mesh.X(),mesh.Y(),np.reshape(sol,(nvy,nvx)))
        plt.title(title)
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        plt.colorbar()
        dim = mesh.dim()
        Nx=mesh.NdomX() ; Ny=mesh.NdomY() ; Nz=mesh.NdomZ()
        tagsCut=np.reshape(tags,(Nz,Ny,Nx))[1:-1,1:-1,1:-1]
        tagsCut=tagsCut.flatten()
        data = {'X':X,'Y':Y,'Z':Z,'Tag':tagsCut,'phi':sol}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        ##graficar primero probar qué está imprimiendo
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]
        
        if f_anlytc:
            X=list(X)
            X=np.array(X)
            exac_phi=f_anlytc(X,Y,Z)
            error = phi-exac_phi
            max_e =np.amax(error)
            distance = (np.dot(error,error))**0.5
        #bool_cond=df['Tag'].str.contains('D', regex=False)
        #dirich_coordTags=np.transpose(df[bool_cond].values)
        #bool_cond=df['Tag'].str.contains('N', regex=False)
        #neum_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------
        if dim ==1:
            plt.figure()
            plt.plot(X,phi,'o')
            plt.title(title)
            plt.xlabel('X')
            plt.ylabel('f(X)')
            plt.xlim(0.,mesh.lx())
            if f_anlytc:
                plt.figure()
                plt.plot(X, 100*error/exac_phi ,'s')
                plt.title('Max Err={:.5f} Dist={:.5f}'.format(max_e,distance))
                plt.xlabel('X')
                plt.ylabel('% Error')
        elif dim ==2:
            #plt.figure()
            plt.scatter(X,Y,c=phi, cmap='gnuplot')
            #plt.contour(mesh.X(),mesh.Y(),np.reshape(sol,(nvy,nvx)))
            #plt.title(title)
            #plt.xlim(0.,mesh.lx())
            #plt.ylim(0.,mesh.ly())
            plt.colorbar()
        elif dim ==3:
            fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.scatter(dirich_coordTags[0], dirich_coordTags[1], dirich_coordTags[2], color='y')
            #ax.scatter(neum_coordTags[0],neum_coordTags[1] ,neum_coordTags[2] , color='b')
            #ax.scatter(source_coordTags[0],source_coordTags[1] , source_coordTags[2], color='m')
            #ax.scatter(inter_coordTags[0],inter_coordTags[1] , inter_coordTags[2], color='k')
        plt.show()   
        
    def energyVel(self,sol,vel,stgd_positions,f_anlytc=None,title='Solution',size=20,maxSpeed=None, colorV='g',marcador="o", barra=True):
        
        mesh=self._mesh
        dim = mesh.dim()
        tags = mesh.tags()
        scalarMesh_X=mesh.X()
        scalarMesh_Y=mesh.Y()
        scalarMesh_Z=mesh.Z()
        coordY,coordZ,coordX = np.meshgrid(scalarMesh_Y,scalarMesh_Z,scalarMesh_X)
        sol = sol.flatten()
        X=coordX.flatten()
        Y=coordY.flatten()
        Z=coordZ.flatten()
        Nx=mesh.NdomX() ; Ny=mesh.NdomY() ; Nz=mesh.NdomZ()
        tagsCut=np.reshape(tags,(Nz,Ny,Nx))[1:-1,1:-1,1:-1]
        tagsCut=tagsCut.flatten()
        #print(X.shape,Y.shape,Z.shape,tagsCut.shape,sol.shape)
        data = {'X':X,'Y':Y,'Z':Z,'Tag':tagsCut,'phi':sol}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I|S', regex=True)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        ##graficar primero probar qué está imprimiendo
        coordTagsSol=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        X = coordTagsSol[0] ; Y = coordTagsSol[1]
        Z = coordTagsSol[2] ; phi = coordTagsSol[4]
        UMesh_X = stgd_positions[0] ; UMesh_Y = scalarMesh_Y ; UMesh_Z = scalarMesh_Z
        U = vel[0]
        if dim>1:
            VMesh_X = scalarMesh_X ; VMesh_Y = stgd_positions[1] ; UMesh_Z = scalarMesh_Z
            V = vel[1]
        if dim>2:
            WMesh_X = scalarMesh_X ; WMesh_Y = scalarMesh_Y  ; UMesh_Z = stgd_positions[2]
            W = vel[2]
            
        #-------------------Interpolación ----------------------------
        #..............¡¡¡SOLO FUNCIONA PARA MALLAS REGULARES 2D!!!.......
        if dim ==2:
            U_inter=0.5*(U[0,:-1,:]+U[0,1:,:])
            V_inter=0.5*(V[0,:,:-1]+V[0,:,1:])
            x_uniform,y_uniform = np.meshgrid( UMesh_X, VMesh_Y) 
        #if dim==3:
        
        if f_anlytc:
            X=list(X)
            X=np.array(X)
            exac_phi=f_anlytc(X,Y,Z)
            error = phi-exac_phi
            max_e =np.amax(error)
            distance = (np.dot(error,error))**0.5
        #bool_cond=df['Tag'].str.contains('D', regex=False)
        #dirich_coordTags=np.transpose(df[bool_cond].values)
        #bool_cond=df['Tag'].str.contains('N', regex=False)
        #neum_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------
        if dim ==1:
            plt.figure()
#            plt.plot(X,phi,'o',markersize=5)
#            plt.title(title)
#            plt.xlabel('X')
#            plt.ylabel('f(X)')
#            plt.xlim(0.,mesh.lx())
#            if f_anlytc:
#                plt.figure()
#                plt.plot(X, 100*error/exac_phi ,'s')
#                plt.title('Max Err={:.5f} Dist={:.5f}'.format(max_e,distance))
#                plt.xlabel('X')
#                plt.ylabel('% Error')
        elif dim ==2:
            plt.figure()
            plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, marker=marcador)
            plt.title(title)
            plt.xlim(0.,mesh.lx())
            plt.ylim(0.,mesh.ly())
            if barra: plt.colorbar()
            speed = np.sqrt(U_inter*U_inter + V_inter*V_inter)
            if maxSpeed: lw = 2.*speed / maxSpeed
            else: lw = 2.*speed / speed.max()
            print("La rapidez promedio es:", np.average(speed))
            plt.streamplot(x_uniform, y_uniform, U_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
            #plt.xlim((0,lx))
            #plt.ylim((0,ly))
            #plt.show()
        elif dim ==3:
            fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.scatter(dirich_coordTags[0], dirich_coordTags[1], dirich_coordTags[2], color='y')
            #ax.scatter(neum_coordTags[0],neum_coordTags[1] ,neum_coordTags[2] , color='b')
            #ax.scatter(source_coordTags[0],source_coordTags[1] , source_coordTags[2], color='m')
            #ax.scatter(inter_coordTags[0],inter_coordTags[1] , inter_coordTags[2], color='k')
        plt.show()
        
    def scalarVelCuts(self,sol,vel,x_cut=None,y_cut=None,z_cut=None,title='Solution',size=20,maxSpeed=None, colorV='g',marcador="o",min_val=None,max_val=None, barra=True):
        
        mesh=self._mesh
        volX=mesh.X() ; volY=mesh.Y() ; volZ=mesh.Z()
        U = vel[0] ; V=vel[1] ; W=vel[2]
        if x_cut==None: x_cut=mesh.lx()/2.
        if y_cut==None: y_cut=mesh.ly()/2.
        if z_cut==None: z_cut=mesh.lz()/2.
        indx_xcut=np.argmax(x_cut<volX)
        indx_ycut=np.argmax(y_cut<volY)
        indx_zcut=np.argmax(z_cut<volZ)
        if (x_cut-volX[indx_xcut-1]) < (volX[indx_xcut] - x_cut): 
            indx_xcut-=1  #Si está mas cerca de X[índice anterior] tomamos el índice anterior
        if (y_cut-volY[indx_ycut-1]) < (volY[indx_ycut] - y_cut): 
            indx_ycut-=1
        if (z_cut-volZ[indx_zcut-1]) < (volZ[indx_zcut] - z_cut): 
            indx_zcut-=1
        
        #-------------------Interpolación ----------------------------
        #..............¡¡¡SOLO FUNCIONA PARA MALLAS REGULARES 2D!!!.......
        #==========================================================
        #           Y CUT
        #============================================================
        phi = sol[:,indx_ycut,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY[indx_ycut],volZ,volX)
        x_uniform,y_uniform = np.meshgrid( volZ, volX) 
        X=X.flatten() ; Z=Z.flatten()
        U_inter=0.5*(U[:,indx_ycut,:-1]+U[:,indx_ycut,1:])
        W_inter=0.5*(W[:-1,indx_ycut,:]+W[1:,indx_ycut,:])
        plt.figure()
        plt.scatter(Z,X,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZX con y={:.3f}".format(volY[indx_ycut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.lx())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + W_inter*W_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
#        print(x_uniform.shape,y_uniform.shape,W_inter.shape, U_inter.shape)
        plt.streamplot(x_uniform, y_uniform, W_inter, U_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()
        #=========================================================
        #           X CUT
        #====================================================
        phi = sol[:,:,indx_xcut]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ,volX[indx_xcut])
        x_uniform,y_uniform = np.meshgrid( volZ, volY) 
        Z=Z.flatten() ; Y=Y.flatten()
        W_inter=0.5*(W[:-1,:,indx_xcut]+W[1:,:,indx_xcut])
        V_inter=0.5*(V[:,:-1,indx_xcut]+V[:,1:,indx_xcut])
        plt.figure()
        plt.scatter(Z,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZY con x={:.3f}".format(volX[indx_xcut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(W_inter*W_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(x_uniform, y_uniform, W_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()
        #==============================================================
        #           Z CUT
        #====================================================
        phi = sol[indx_zcut,:,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ[indx_zcut],volX)
        x_uniform,y_uniform = np.meshgrid( volX, volY) 
        X=X.flatten() ; Y=Y.flatten()
        U_inter=0.5*(U[indx_zcut,:,:-1]+U[indx_zcut,:,1:])
        V_inter=0.5*(V[indx_zcut,:-1,:]+V[indx_zcut,1:,:])
        plt.figure()
        plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano XY con z={:.3f}".format(volZ[indx_zcut]) )
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(x_uniform, y_uniform, U_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()
        
    def scalarVelCutsColor(self,sol,vel,x_cut=None,y_cut=None,z_cut=None,title='Solution',size=20,maxSpeed=None, colorV='g',marcador="o",min_val=None,max_val=None, barra=True, barraX=False,barraY=False,barraZ=False):
        
        mesh=self._mesh
        volX=mesh.X() ; volY=mesh.Y() ; volZ=mesh.Z()
        U = vel[0] ; V=vel[1] ; W=vel[2]
        if x_cut==None: x_cut=mesh.lx()/2.
        if y_cut==None: y_cut=mesh.ly()/2.
        if z_cut==None: z_cut=mesh.lz()/2.
        indx_xcut=np.argmax(x_cut<volX)
        indx_ycut=np.argmax(y_cut<volY)
        indx_zcut=np.argmax(z_cut<volZ)
        if (x_cut-volX[indx_xcut-1]) < (volX[indx_xcut] - x_cut): 
            indx_xcut-=1  #Si está mas cerca de X[índice anterior] tomamos el índice anterior
        if (y_cut-volY[indx_ycut-1]) < (volY[indx_ycut] - y_cut): 
            indx_ycut-=1
        if (z_cut-volZ[indx_zcut-1]) < (volZ[indx_zcut] - z_cut): 
            indx_zcut-=1
        
        #-------------------Interpolación ----------------------------
        #..............¡¡¡SOLO FUNCIONA PARA MALLAS REGULARES 2D!!!.......
        #==========================================================
        #           Y CUT
        #============================================================
        phi = sol[:,indx_ycut,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY[indx_ycut],volZ,volX)
        x_uniform,y_uniform = np.meshgrid( volZ, volX) 
        X=X.flatten() ; Z=Z.flatten()
        U_inter=0.5*(U[:,indx_ycut,:-1]+U[:,indx_ycut,1:])
        W_inter=0.5*(W[:-1,indx_ycut,:]+W[1:,indx_ycut,:])
        plt.figure()
        plt.scatter(Z,X,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZX con y={:.3f}".format(volY[indx_ycut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.lx())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + W_inter*W_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
#        print(x_uniform.shape,y_uniform.shape,W_inter.shape, U_inter.shape)
        plt.streamplot(x_uniform, y_uniform, W_inter, U_inter, density=0.6, color=colorV, linewidth=lw)
        vis.color2D((mesh.dominioZ(),mesh.dominioX()),phi,'inferno',barra=barraY,turn=True,vm=min_val,vM=max_val)
        plt.show()
        #=========================================================
        #           X CUT
        #====================================================
        phi = sol[:,:,indx_xcut]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ,volX[indx_xcut])
        x_uniform,y_uniform = np.meshgrid( volZ, volY) 
        Z=Z.flatten() ; Y=Y.flatten()
        W_inter=0.5*(W[:-1,:,indx_xcut]+W[1:,:,indx_xcut])
        V_inter=0.5*(V[:,:-1,indx_xcut]+V[:,1:,indx_xcut])
        plt.figure()
        plt.scatter(Z,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZY con x={:.3f}".format(volX[indx_xcut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(W_inter*W_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(x_uniform, y_uniform, W_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        vis.color2D((mesh.dominioZ(),mesh.dominioY()),phi,'inferno', barra=barraX, turn=True,vm=min_val,vM=max_val)
        plt.show()
        #==============================================================
        #           Z CUT
        #====================================================
        phi = sol[indx_zcut,:,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ[indx_zcut],volX)
        x_uniform,y_uniform = np.meshgrid( volX, volY) 
        X=X.flatten() ; Y=Y.flatten()
        U_inter=0.5*(U[indx_zcut,:,:-1]+U[indx_zcut,:,1:])
        V_inter=0.5*(V[indx_zcut,:-1,:]+V[indx_zcut,1:,:])
        plt.figure()
        plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano XY con z={:.3f}".format(volZ[indx_zcut]) )
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(x_uniform, y_uniform, U_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        vis.color2D((mesh.dominioX(),mesh.dominioY()),phi,'inferno',barra=barraZ,vm=min_val,vM=max_val)
        plt.show()

    def scalarVelCutswError(self,sol,vel,x_cut=None,y_cut=None,z_cut=None,title='Solution',size=20,maxSpeed=None, colorV='g',marcador="o",min_val=None,max_val=None, barra=True):
        
        mesh=self._mesh
        volX=mesh.X() ; volY=mesh.Y() ; volZ=mesh.Z()
        U = vel[0] ; V=vel[1] ; W=vel[2]
        if x_cut==None: x_cut=mesh.lx()/2.
        if y_cut==None: y_cut=mesh.ly()/2.
        if z_cut==None: z_cut=mesh.lz()/2.
        indx_xcut=np.argmax(x_cut<volX)
        indx_ycut=np.argmax(y_cut<volY)
        indx_zcut=np.argmax(z_cut<volZ)
        if (x_cut-volX[indx_xcut-1]) < (volX[indx_xcut] - x_cut): 
            indx_xcut-=1  #Si está mas cerca de X[índice anterior] tomamos el índice anterior
        if (y_cut-volY[indx_ycut-1]) < (volY[indx_ycut] - y_cut): 
            indx_ycut-=1
        if (z_cut-volZ[indx_zcut-1]) < (volZ[indx_zcut] - z_cut): 
            indx_zcut-=1
        
        #-------------------Interpolación ----------------------------
        #..............¡¡¡SOLO FUNCIONA PARA MALLAS REGULARES 2D!!!.......
        #==========================================================
        #           Y CUT
        #============================================================
        phi = sol[:,indx_ycut,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY[indx_ycut],volZ,volX)
        x_uniform,z_uniform = np.meshgrid( volX, volZ) 
        X=X.flatten() ; Z=Z.flatten()
        U_inter=0.5*(U[:,indx_ycut,:-1]+U[:,indx_ycut,1:])
        W_inter=0.5*(W[:-1,indx_ycut,:]+W[1:,indx_ycut,:])
        plt.figure()
        plt.scatter(Z,X,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZX con y={:.3f}".format(volY[indx_ycut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.lx())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + W_inter*W_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        print(x_uniform.shape,z_uniform.shape,W_inter.shape, U_inter.shape)
        plt.streamplot(z_uniform, x_uniform, W_inter, U_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()
        #=========================================================
        #           X CUT
        #====================================================
        phi = sol[:,:,indx_xcut]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ,volX[indx_xcut])
        z_uniform,y_uniform = np.meshgrid( volZ, volY) 
        Z=Z.flatten() ; Y=Y.flatten()
        W_inter=0.5*(W[:-1,:,indx_xcut]+W[1:,:,indx_xcut])
        V_inter=0.5*(V[:,:-1,indx_xcut]+V[:,1:,indx_xcut])
        plt.figure()
        plt.scatter(Z,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano ZY con x={:.3f}".format(volX[indx_xcut]) )
        plt.xlim(0.,mesh.lz())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(W_inter*W_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(z_uniform, y_uniform, W_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()
        #==============================================================
        #           Z CUT
        #====================================================
        phi = sol[indx_zcut,:,:]
        phi = phi.flatten()
        Y,Z,X = np.meshgrid(volY,volZ[indx_zcut],volX)
        x_uniform,y_uniform = np.meshgrid( volX, volY) 
        X=X.flatten() ; Y=Y.flatten()
        U_inter=0.5*(U[indx_zcut,:,:-1]+U[indx_zcut,:,1:])
        V_inter=0.5*(V[indx_zcut,:-1,:]+V[indx_zcut,1:,:])
        plt.figure()
        plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, vmin=min_val,vmax=max_val, marker=marcador)
        plt.title( title+" Plano XY con z={:.3f}".format(volZ[indx_zcut]) )
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        if barra: plt.colorbar()
        speed = np.sqrt(U_inter*U_inter + V_inter*V_inter)
        if maxSpeed: lw = 2.*speed / maxSpeed
        else: lw = 2.*speed / speed.max()
        print("La rapidez promedio es:", np.average(speed))
        plt.streamplot(x_uniform, y_uniform, U_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
        plt.show()


    def scalarVel(self,sol,vel,stgd_positions,f_anlytc=None,title='Solution',size=20,maxSpeed=None, colorV='g',marcador="o", barra=True):
        
        mesh=self._mesh
        dim = mesh.dim()
        tags = mesh.tags()
        scalarMesh_X=mesh.X()
        scalarMesh_Y=mesh.Y()
        scalarMesh_Z=mesh.Z()
        coordY,coordZ,coordX = np.meshgrid(scalarMesh_Y,scalarMesh_Z,scalarMesh_X)
        sol = sol.flatten()
        X=coordX.flatten()
        Y=coordY.flatten()
        Z=coordZ.flatten()
        Nx=mesh.NdomX() ; Ny=mesh.NdomY() ; Nz=mesh.NdomZ()
        tagsCut=np.reshape(tags,(Nz,Ny,Nx))[1:-1,1:-1,1:-1]
        tagsCut=tagsCut.flatten()
        phi = sol
        UMesh_X = stgd_positions[0] ; UMesh_Y = scalarMesh_Y ; UMesh_Z = scalarMesh_Z
        U = vel[0]
        if dim>1:
            VMesh_X = scalarMesh_X ; VMesh_Y = stgd_positions[1] ; UMesh_Z = scalarMesh_Z
            V = vel[1]
        if dim>2:
            WMesh_X = scalarMesh_X ; WMesh_Y = scalarMesh_Y  ; UMesh_Z = stgd_positions[2]
            W = vel[2]
            
        #-------------------Interpolación ----------------------------
        #..............¡¡¡SOLO FUNCIONA PARA MALLAS REGULARES 2D!!!.......
        if dim ==2:
            U_inter=0.5*(U[0,:-1,:]+U[0,1:,:])
            V_inter=0.5*(V[0,:,:-1]+V[0,:,1:])
            x_uniform,y_uniform = np.meshgrid( UMesh_X, VMesh_Y) 
        #if dim==3:
        
        if f_anlytc:
            X=list(X)
            X=np.array(X)
            exac_phi=f_anlytc(X,Y,Z)
            error = phi-exac_phi
            max_e =np.amax(error)
            distance = (np.dot(error,error))**0.5
        #bool_cond=df['Tag'].str.contains('D', regex=False)
        #dirich_coordTags=np.transpose(df[bool_cond].values)
        #bool_cond=df['Tag'].str.contains('N', regex=False)
        #neum_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------
        if dim ==1:
            plt.figure()
#            plt.plot(X,phi,'o',markersize=5)
#            plt.title(title)
#            plt.xlabel('X')
#            plt.ylabel('f(X)')
#            plt.xlim(0.,mesh.lx())
#            if f_anlytc:
#                plt.figure()
#                plt.plot(X, 100*error/exac_phi ,'s')
#                plt.title('Max Err={:.5f} Dist={:.5f}'.format(max_e,distance))
#                plt.xlabel('X')
#                plt.ylabel('% Error')
        elif dim ==2:
            plt.figure()
            plt.scatter(X,Y,c=phi, cmap='gnuplot', s=size, marker=marcador)
            plt.title(title)
            plt.xlim(0.,mesh.lx())
            plt.ylim(0.,mesh.ly())
            if barra: plt.colorbar()
            speed = np.sqrt(U_inter*U_inter + V_inter*V_inter)
            if maxSpeed: lw = 2.*speed / maxSpeed
            else: lw = 2.*speed / speed.max()
            print("La rapidez promedio es:", np.average(speed))
            plt.streamplot(x_uniform, y_uniform, U_inter, V_inter, density=0.6, color=colorV, linewidth=lw)
            #plt.xlim((0,lx))
            #plt.ylim((0,ly))
            #plt.show()
        elif dim ==3:
            fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            #ax.scatter(dirich_coordTags[0], dirich_coordTags[1], dirich_coordTags[2], color='y')
            #ax.scatter(neum_coordTags[0],neum_coordTags[1] ,neum_coordTags[2] , color='b')
            #ax.scatter(source_coordTags[0],source_coordTags[1] , source_coordTags[2], color='m')
            #ax.scatter(inter_coordTags[0],inter_coordTags[1] , inter_coordTags[2], color='k')
        plt.show()
        
    def Contour(self,sol, tieneNiv=False ,niveles=None ,title='Solution'):
        
        mesh=self._mesh
        tags = mesh.tags()
        volX=mesh.X()
        volY=mesh.Y()
        volZ=mesh.Z()
        nvx =mesh.nvx() ; nvy = mesh.nvy()
        coordY,coordZ,coordX = np.meshgrid(volY,volZ,volX)
        sol=sol.flatten()
        X=coordX.flatten()
        Y=coordY.flatten()
        Z=coordZ.flatten()
        plt.figure()
        if not(tieneNiv): CS=plt.contour(mesh.X(),mesh.Y(),np.reshape(sol,(nvy,nvx)))
        else:  CS=plt.contour(mesh.X(),mesh.Y() , np.reshape(sol,(nvy,nvx)), niveles )
        plt.clabel(CS)
        plt.title(title)
        plt.xlim(0.,mesh.lx())
        plt.ylim(0.,mesh.ly())
        #plt.colorbar()
        plt.show()