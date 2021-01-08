# -*- coding: utf-8 -*-
"""
Created on Fri May 24 14:23:43 2019

@author: luiggi
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
import matplotlib.pyplot as plt

def colormap2D(mesh, sol, colormap, vis='con', complete = False):
    X = mesh.dominioX()
    Y = mesh.dominioY()
    nvx = mesh.nvx()
    nvy = mesh.nvy()

    sol.shape = (nvy,nvx)
    
    if complete:
        u = np.zeros((nvx+2,nvy+2))
        u[0,1:-1] = sol[0,:]
        u[-1,1:-1] = sol[-1,:]
        u[:,0] = 100
        u[:,-1] = 500
        u[1:nvx+1,1:nvy+1] = sol
        xg, yg = np.meshgrid(X,Y)
    else:
        u = sol
        xg, yg = np.meshgrid(X[1:-1],Y[1:-1])
    
    if vis == 'con':
        pl.contourf(xg, yg, u, 50, alpha=.95, cmap=colormap)
        pl.colorbar()
#    C = pl.contour(xg, yg, u, 10, colors='black', alpha=0.0, linewidth=.1)
#    pl.clabel(C, inline=1, fontsize=10)
    elif vis == 'surf':
        fig = pl.figure()
        ax = Axes3D(fig)
        ax.plot_surface(xg, yg, u, rstride=2, cstride=2, alpha=.95, cmap=colormap)

    pl.show()
    
    sol.shape = sol.size
    
def color2D(unique_postns, sol, colormap, vis='con', barra=False,complete = False,vm=None,vM=None,turn=False):
    X , Y = unique_postns
    nvx = len(X)-2
    nvy = len(Y)-2

    sol.shape = (nvy,nvx)
    
    if complete:
        u = np.zeros((nvx+2,nvy+2))
        u[0,1:-1] = sol[0,:]
        u[-1,1:-1] = sol[-1,:]
        u[:,0] = 100
        u[:,-1] = 500
        u[1:nvx+1,1:nvy+1] = sol
        xg, yg = np.meshgrid(X,Y)
    else:
        if turn:
            u=np.copy(sol).transpose()
        else:
            u = sol
        xg, yg = np.meshgrid(X[1:-1],Y[1:-1])
    
    if vis == 'con':
        pl.contourf(xg, yg, u, 50, alpha=.95, vmin=vm , vmax=vM ,cmap=colormap)
        if barra: pl.colorbar()
#    C = pl.contour(xg, yg, u, 10, colors='black', alpha=0.0, linewidth=.1)
#    pl.clabel(C, inline=1, fontsize=10)
    elif vis == 'surf':
        fig = pl.figure()
        ax = Axes3D(fig)
        ax.plot_surface(xg, yg, u, rstride=2, cstride=2, alpha=.95, cmap=colormap)

    pl.show()
    
    sol.shape = sol.size

def colormap2D_adv(mesh, sol, colormap, vis='con', complete = False):
    X = mesh.dominioX()
    Y = mesh.dominioY()

    print('SHAPE', sol.shape)
    xg, yg = np.meshgrid(X[1:-1],Y[1:-1])
    
    pl.contourf(xg, yg, sol[0], 50, alpha=.95, cmap=colormap)
#    C = pl.contour(xg, yg, u, 10, colors='black', alpha=0.0, linewidth=.1)
#    pl.clabel(C, inline=1, fontsize=10)

    pl.show()
    