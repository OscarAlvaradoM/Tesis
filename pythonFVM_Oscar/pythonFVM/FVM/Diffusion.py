#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 21:37:18 2019

@author: jose
"""

import numpy as np
#from Coefficients import Coefficients

#def funcionConst(x,y,z):
#    return 10.

class Diffusion():
    """
    Class that has the methods for calculating the modifications of the aP, aE, aW, ...
    coefficients due to the Diffusion term in the transport equation.
    """
    
    def __init__(self, mesh, funcionGamma):
        """ 
        Saves the mesh and a function as attributes of the object.  
        
        funcionGamma: function that can be evaluated to obtain the value of the difussion
        coeficient, Gamma, in a particular position, that is Gamma = f(x,y,z). If 
        funcionGamma is an int or float it is assumed that the difussion coefficient 
        is constant. (float/int/function)
        """
        
        self._mesh = mesh
        self._gammaconstante = 0.
        self._Gamma = None
        
        if isinstance(funcionGamma,(int,float)):
            self._gammaconstante = funcionGamma
            self._Gamma = self.funcionConst
        else:
            self._Gamma = funcionGamma
          
    
    def funcionConst(self, x, y, z):
        return self._gammaconstante
    
    def DE(self):
        """
        Gives a 3D numpy array that should be used to correct the aE coefficient 
        caused by the diffusion effect.
        """
        
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dxe = mesh.dxe() # "Separations between the corresponding node and it's east neighbor." 
        Gamma_e = self._Gamma(mesh.xe(x, dxe), y , z) # xe() -> x coordinate of the east border in the volume.
        DE = Gamma_e * mesh.areasX() / dxe # areasX() -> "the areas of the faces in the X direction."
        return DE

    def DW(self):
        """
        Gives a 3D numpy array that should be used to correct the aW coefficient
        caused by the diffusion effect.
        """        
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dxw = mesh.dxw() # "Separations between the corresponding node and it's west neighbor."
        Gamma_w = self._Gamma(mesh.xw(x, dxw), y , z) # xw() -> "x coordinate of the west border in the volume."
        DW = Gamma_w * mesh.areasX() / dxw # areasX() -> "the areas of the faces in the X direction."
        return DW

    def DN(self):
        """
        Gives a 3D numpy array that should be used to correct the aN coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dyn = mesh.dyn() # "Separations between the corresponding node and it's north neighbor."
        Gamma_n = self._Gamma(x, mesh.yn(y, dyn), z) # yn() -> "y coordinate of the north border in the volume."
        DN = Gamma_n * mesh.areasY() / dyn # areasY() -> "the areas of the faces in the Y direction."
        return DN

    def DS(self):
        """
        Gives a 3D numpy array that should be used to correct the aS coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dys = mesh.dys() # "Separations between the corresponding node and it's south neighbor."
        Gamma_s = self._Gamma(x, mesh.ys(y, dys), z) # ys() -> "y coordinate of the south border in the volume."
        DS = Gamma_s * mesh.areasY() / dys # areasY() -> "the areas of the faces in the Y direction."
        return DS

    def DT(self):
        """
        Gives a 3D numpy array that should be used to correct the aT coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dzt = mesh.dzt() # "Separations between the corresponding node and it's top neighbor."
        Gamma_t = self._Gamma(x, y, mesh.zt(z, dzt)) # zt() -> "z coordinate of the top border in the volume."
        DT = Gamma_t * mesh.areasZ() / dzt # areasZ() -> "the areas of the faces in the Z direction."
        return DT

    def DB(self):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        (y, z, x) = np.meshgrid(mesh.Y(), mesh.Z(), mesh.X())
        dzb = mesh.dzb() # "Separations between the corresponding node and it's bottom neighbor."
        Gamma_b = self._Gamma(x, y, mesh.zb(z, dzb)) # zt() -> "z coordinate of the top border in the volume."
        DB = Gamma_b * mesh.areasZ() / dzb # areasZ() -> "the areas of the faces in the Z direction."
        return DB
