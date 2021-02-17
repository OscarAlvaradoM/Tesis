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
    
    def __init__(self, mesh, gamma):
        """ 
        Saves the mesh and a function as attributes of the object.  
        
        gamma: function that can be evaluated to obtain the value of the difussion
        coeficient, Gamma, in a particular position, that is Gamma = f(x,y,z). If 
        funcionGamma is an int or float it is assumed that the difussion coefficient 
        is constant. (float/int/function)
        """
        
        self._mesh = mesh
        self._gammaconstante = 0.
        self._gamma = None
        
        if isinstance(gamma, (int,float)):
            self._gammaconstante = gamma
            self._gamma = self.funcionConst
        # Se asume que es una función en este caso
        else:
            self._gamma = gamma
          
    
    def funcionConst(self, x, y, z):
        return self._gammaconstante
    
    def east_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aE coefficient 
        caused by the diffusion effect.
        """
        
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_x = mesh.faces[0]
        δ_x = np.array(mesh.deltas[0] + (faces_x[-1] - mesh.coords[0][-1],))
        gamma_e = self._gamma(faces_x[1:], y , z) # Usando todas las de la derecha (east)
        areas = mesh.areas_x()[1:]
        #areas[-1,:,:] = 0
        east_d = gamma_e * areas / δ_x[:,None,None]
        return east_d

    def west_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aW coefficient
        caused by the diffusion effect.
        """        
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_x = mesh.faces[0] 
        δ_x = np.array((mesh.coords[0][0],) + mesh.deltas[0])
        gamma_w = self._gamma(faces_x[:-1], y , z) # Usando todas las de la derecha (east)
        areas = mesh.areas_x()[:-1]
        #areas[0,:,:] = 0
        west_d = gamma_w * areas / δ_x[:,None,None]
        return west_d

    def north_difussion(self):
        """
        Gives a 3D numpy array that should be used to correct the aN coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_y = mesh.faces[1]
        δ_y = np.array(mesh.deltas[1] + (faces_y[-1] - mesh.coords[1][-1],))
        gamma_n = self._gamma(x, faces_y[1:] , z) # Usando todas las de la derecha (east)
        areas = mesh.areas_y()[1:]
        #areas[-1,:,:] = 0
        north_d = gamma_n * areas / δ_y[None,:,None]
        return north_d

    def south_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aS coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_y = mesh.faces[1] 
        δ_y = np.array((mesh.coords[1][0],) + mesh.deltas[1])
        gamma_s = self._gamma(x, faces_y[:-1] , z) # Usando todas las de la derecha (east)
        areas = mesh.areas_y()[:-1]
        #areas[0,:,:] = 0
        south_d = gamma_s * areas / δ_y[None,:,None]
        return south_d

    def top_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aT coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_z = mesh.faces[2]
        δ_z = np.array(mesh.deltas[2] + (1,))
        gamma_t = self._gamma(x, y, faces_z[1:]) # Usando todas las de la derecha (east)
        areas = mesh.areas_z()[1:]
        #areas[-1,:,:] = 0
        top_d = gamma_t * areas / δ_y[None,None,:]
        return top_d

    def bottom_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        faces_z = mesh.faces[2]
        δ_z = np.array((1,) + mesh.deltas[2])
        gamma_b = self._gamma(x, y, faces_z[:-1]) # Usando todas las de la derecha (east)
        areas = mesh.areas_z()[:-1]
        #areas[0,:,:] = 0
        bottom_d = gamma_b * areas / δ_z[None,None,:]
        return bottom_d
