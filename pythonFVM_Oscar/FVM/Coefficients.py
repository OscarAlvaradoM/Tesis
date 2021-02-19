#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 21:14:10 2019

@author: jose
"""

import numpy as np
#from numpy import maximum as maxi
from Diffusion import Diffusion
#import scipy.sparse as sp

class Coefficients():
    """
    Class that defines the coefficients that every node needs to form it's
    discrete equation that represent the diferential equation. All the coefficients
    needed for a scalar quantity (also called variable or unkown) are contained
    in one object (instance of this class).
    """    

    def __init__(self, mesh):
        self._dim = mesh.dim
        self._mesh = mesh
        self._volumes = mesh.volumes
        self._Su = np.zeros(self._volumes)  ;  self._Sp = None
        self._aN = None ; self._aS = None ; self._aT = None ; self._aB = None
        self.init_coefs()
        
    def init_coefs(self):
        """
        Método que inicializa los coeficientes en un arreglo tridimensional
        """
        vols = self._volumes
        self._aP = np.zeros(vols)
        dim = self._dim
        self._aE = np.zeros(vols)
        self._aW = np.zeros(vols)
        self._Sp = np.zeros(vols)
        self._Su = np.zeros(vols)
        if dim > 1:
            self._aN = np.zeros(vols)
            self._aS = np.zeros(vols)
        if dim == 3:
            self._aT = np.zeros(vols)
            self._aB = np.zeros(vols)
        
    def set_diffusion(self, gamma):
        """This method makes an instance of the Diffusion class and uses it to
        update the aP,aE,aW, ... atributes.
        
        Gamma: The diffusion coeffcient of the transport equation, if the diffusion
               coefficient is constant then a float can be given but in other cases
               a function of the coordinates, denoted as Gamma=f(x,y,z) should be 
               used. (float or function)"""
        
        dim = self._dim
        malla = self._mesh
        diffusion = Diffusion(malla, gamma)
        east_diff, sp_e, su_e = diffusion.east_diffusion()
        west_diff, sp_w, su_w = diffusion.west_diffusion()
        self._aE -= east_diff
        self._aW -= west_diff
        self._Sp -= sp_e + sp_w
        self._Su += su_e + su_w
        self._aP +=  self._aE + self._aW
        if dim > 1:
            north_diff, sp_n, su_n = diffusion.north_diffusion()
            south_diff, sp_s, su_s = diffusion.south_diffusion()
            self._aN -= north_diff
            self._aS -= south_diff
            self._Sp -= sp_n + sp_s
            self._Su += su_n + su_s
            self._aP += self._aN + self._aS
        if dim == 3:
            top_diff, sp_t, su_t = diffusion.top_diffussion()
            bottom_diff, sp_b, su_b = diffusion.bottom_diffusion()
            self._aT -= top_diff
            self._aB -= bottom_diff
            self._Sp -= sp_t + sp_b
            self._Su += su_t + su_b
            self._aP += top_diff + bottom_diff
        self._aP += self._Sp
            
    def get_Su(self):
        return self._Su
    
    def get_Sp(self):
        return self._Sp
    
    def get_aP(self):
        return self._aP
    
    def get_aE(self):
        return self._aE
    
    def get_aW(self):
        return self._aW

    def get_aN(self):
        return self._aN

    def get_aS(self):
        return self._aS

    def get_aT(self):
        return self._aT

    def get_aB(self):
        return self._aB