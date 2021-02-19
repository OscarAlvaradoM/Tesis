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
    
    def east_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aE coefficient 
        caused by the diffusion effect.
        """
        mesh = self._mesh
        faces_x = mesh.faces[0]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=0, orientation="e")
        δ_x = δ_d_x_grid
        gamma_e = self._gamma(faces_x[1:], y , z) # Usando todas las de la derecha (east)
        areas = mesh.areas_x()
        source = areas.copy()
        areas[-1,:,:] = 0
        source[:-1,:,:] = 0
        east_d = gamma_e * areas / δ_x
        # Aquí obtenemos Sp
        sp = gamma_e * source / δ_x
        condicion = malla.get_mask_boundaries_Sp(direction="E")
        sp[-1:,:,:] = sp[-1:,:,:]*np.array(condicion).reshape(sp[-1:,:,:].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="E", gamma = self._gamma)
        su = source
        su[-1:,:,:] = su[-1:,:,:]*np.array(conds_su).reshape(su[-1:,:,:].shape)
        div = δ_x[-1:,:,:]*np.array(condicion).reshape(su[-1:,:,:].shape)
        div = np.where(div == 0., 1, div)
        su[-1:,:,:] = su[-1:,:,:]/div
        return east_d, sp, su

    def west_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aW coefficient
        caused by the diffusion effect.
        """        
        mesh = self._mesh
        faces_x = mesh.faces[0]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=0, orientation="w")
        δ_x = δ_d_x_grid
        gamma_w = self._gamma(faces_x[:-1], y , z) # Usando todas las de la izquierda (west)
        areas = mesh.areas_x()
        source = areas.copy()
        areas[0,:,:] = 0
        source[1:,:,:] = 0
        west_d = gamma_w * areas / δ_x
        # Aquí obtenemos Sp
        sp = gamma_w * source / δ_x
        condicion = malla.get_mask_boundaries_Sp(direction="W")
        sp[:1,:,:] = sp[:1,:,:]*np.array(condicion).reshape(sp[:1,:,:].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="W", gamma = self._gamma)
        su = source
        su[:1,:,:] = su[:1,:,:]*np.array(conds_su).reshape(su[:1,:,:].shape)
        div = δ_x[:1,:,:]*np.array(condicion).reshape(su[:1,:,:].shape)
        div = np.where(div == 0., 1, div)
        su[:1,:,:] = su[:1,:,:]/div
        return west_d, sp, su

    def north_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aN coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        faces_y = mesh.faces[1]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=1, orientation="n")
        δ_y = δ_d_y_grid
        gamma_n = self._gamma(x, faces_y[1:] , z) # Usando todas las de arriba (north)
        areas = mesh.areas_y()
        source = areas.copy()
        areas[:,-1,:] = 0
        source[:,:-1,:] = 0
        north_d = gamma_n * areas / δ_y
        # Aquí obtenemos Sp
        sp = gamma_n * source / δ_y
        condicion = malla.get_mask_boundaries_Sp(direction="N")
        sp[:,-1:,:] = sp[:,-1:,:]*np.array(condicion).reshape(sp[:,-1:,:].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="N", gamma = self._gamma)
        su = source
        su[:,-1:,:] = su[:,-1:,:]*np.array(conds_su).reshape(su[:,-1:,:].shape)
        div = δ_y[:,-1:,:]*np.array(condicion).reshape(su[:,-1:,:].shape)
        div = np.where(div == 0., 1, div)
        su[:,-1:,:] = su[:,-1:,:]/div
        return north_d, sp, su

    def south_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aS coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        faces_y = mesh.faces[1]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=1, orientation="s")
        δ_y = δ_d_y_grid
        gamma_s = self._gamma(x, faces_y[:-1] , z) # Usando todas las de abajo (south)
        areas = mesh.areas_y()
        source = areas.copy()
        areas[:,0,:] = 0
        source[:,1:,:] = 0
        south_d = gamma_s * areas / δ_y
        # Aquí obtenemos Sp
        sp = gamma_s * source / δ_y
        condicion = malla.get_mask_boundaries_Sp(direction="S")
        sp[:,:1,:] = sp[:,:1,:]*np.array(condicion).reshape(sp[:,:1,:].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="S", gamma = self._gamma)
        su = source
        su[:,:1,:] = su[:,:1,:]*np.array(conds_su).reshape(su[:,:1,:].shape)
        div = δ_y[:,:1,:]*np.array(condicion).reshape(su[:,:1,:].shape)
        div = np.where(div == 0., 1, div)
        su[:,:1,:] = su[:,:1,:]/div
        return south_d, sp, su

    def top_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aT coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        faces_z = mesh.faces[2]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=2, orientation="t")
        δ_z = δ_d_z_grid
        gamma_t = self._gamma(x, y , faces_z[1:]) # Usando todas las superiores (top)
        areas = mesh.areas_z()
        source = areas.copy()
        areas[:,:,-1] = 0
        source[:,:,:-1] = 0
        top_d = gamma_t * areas / δ_z
        # Aquí obtenemos Sp
        sp = gamma_t * source / δ_z
        condicion = malla.get_mask_boundaries_Sp(direction="T")
        sp[:,:,-1:] = sp[:,:,-1:]*np.array(condicion).reshape(sp[:,:,-1:].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="T", gamma = self._gamma)
        su = source
        su[:,:,-1:] = su[:,:,-1:]*np.array(conds_su).reshape(su[:,:,-1:].shape)
        print("Su_t:", su, su.shape)
        print("δ_z_t", δ_z, δ_z.shape)
        div = δ_z[:,:,-1:]*np.array(condicion).reshape(su[:,:,-1:].shape)
        div = np.where(div == 0., 1, div)
        su[:,:,-1:] = su[:,:,-1:]/div
        return top_d, sp, su

    def bottom_diffusion(self):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the diffusion effect.
        """
        mesh = self._mesh
        faces_z = mesh.faces[2]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])
        δ_d_x_grid, δ_d_y_grid, δ_d_z_grid = mesh.get_grid_deltas_dominios(axis=2, orientation="b")
        δ_z = δ_d_z_grid
        gamma_b = self._gamma(x, y, faces_z[:-1]) # Usando todas las inferiores (bottom)
        areas = mesh.areas_z()
        source = areas.copy()
        areas[:,:,0] = 0
        source[:,:,1:] = 0
        bottom_d = gamma_b * areas / δ_z
        # Aquí obtenemos Sp
        sp = gamma_b * source / δ_z
        condicion = malla.get_mask_boundaries_Sp(direction="B")
        sp[:,:,:1] = sp[:,:,:1]*np.array(condicion).reshape(sp[:,:,:1].shape)
        
        # Aquí obtenemos Su
        conds_su = malla.get_mask_boundaries_Su(direction="B", gamma = self._gamma)
        su = source
        su[:,:,:1] = su[:,:,:1]*np.array(conds_su).reshape(su[:,:,:1].shape)
        print("Su_b:", su, su.shape)
        print("δ_z_b", δ_z, δ_z.shape)
        div = δ_z[:,:,:1]*np.array(condicion).reshape(su[:,:,:1].shape)
        div = np.where(div == 0., 1, div)
        su[:,:,:1] = su[:,:,:1]/div
        return bottom_d, sp, su