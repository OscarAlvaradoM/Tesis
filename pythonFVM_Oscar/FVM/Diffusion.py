#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 21:37:18 2019

@author: jose
"""

import numpy as np
#from Coefficients import Coefficients

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
    
    def get_diffusion_coef(self, direction):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the diffusion effect.
        """
        idx, idx_1 = 0,0
        if direction in ["E", "N", "T"]: idx_1=-1
        if direction in ["N", "S"]: idx=1
        elif direction in ["T", "B"]: idx = 2
        mesh = self._mesh
        faces = mesh.faces[idx]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        
        # Para obtener los coeficientes
        δ_d_grid = mesh.get_grid_deltas_dominios(axis=idx, orientation=direction)
        δ = δ_d_grid[idx]
        λ_1 = lambda x, d: x[1:] if d in ["E","N","T"] else x[:-1]
        coord = [var if i!=idx else λ_1(faces, direction) for (i, var) in enumerate([x,y,z])]
        gamma = self._gamma(coord[0], coord[1], coord[2])
        areas = mesh.get_area(direction=idx)
        source = areas.copy()
        coord_2 = [slice(None,None) if var!=idx else idx_1 for var in range(3)]
        areas[coord_2[0],coord_2[1],coord_2[2]] = 0
        λ_2 = lambda d: slice(None,-1) if d in ["E","N","T"] else slice(1,None)
        coord_3 = [slice(None,None) if var!=idx else λ_2(direction) for var in range(3)]
        source[coord_3[0],coord_3[1],coord_3[2]] = 0
        diff = gamma * areas / δ

        # Aquí obtenemos Sp
        sp = gamma * source / δ
        condicion = mesh.get_mask_boundaries_Sp(direction)
        λ_3 = lambda d: slice(-1,None) if d in ["E","N","T"] else slice(None,1)
        coord_4 = [slice(None,None) if var!=idx else λ_3(direction) for var in range(3)]
        tmp_sp = sp[coord_4[0],coord_4[1],coord_4[2]]
        sp[coord_4[0],coord_4[1],coord_4[2]] = tmp_sp*np.array(condicion).reshape(tmp_sp.shape)
        
        # Aquí obtenemos Su
        conds_su = mesh.get_mask_boundaries_Su(direction, gamma = self._gamma)
        su = source
        tmp_su = su[coord_4[0],coord_4[1],coord_4[2]]
        su[coord_4[0],coord_4[1],coord_4[2]] = tmp_su*np.array(conds_su).reshape(tmp_su.shape)
        tmp_su = su[coord_4[0],coord_4[1],coord_4[2]]
        div = δ[coord_4[0],coord_4[1],coord_4[2]]*np.array(condicion).reshape(tmp_su.shape)
        div = np.where(div == 0., 1, div)
        su[coord_4[0],coord_4[1],coord_4[2]] = tmp_su/div
        return diff, sp, su
    
    
    def east_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aE coefficient 
        caused by the diffusion effect.
        """
        east_d, sp, su = self.get_diffusion_coef("E")
        return east_d, sp, su

    def west_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aW coefficient
        caused by the diffusion effect.
        """        
        west_d, sp, su = self.get_diffusion_coef("W")
        return west_d, sp, su

    def north_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aN coefficient
        caused by the diffusion effect.
        """
        north_d, sp, su = self.get_diffusion_coef("N")
        return north_d, sp, su

    def south_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aS coefficient
        caused by the diffusion effect.
        """
        south_d, sp, su = self.get_diffusion_coef("S")
        return south_d, sp, su

    def top_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aT coefficient
        caused by the diffusion effect.
        """
        top_d, sp, su = self.get_diffusion_coef("T")
        return top_d, sp, su

    def bottom_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the diffusion effect.
        """
        bottom_d, sp, su = self.get_diffusion_coef("B")
        return bottom_d, sp, su
