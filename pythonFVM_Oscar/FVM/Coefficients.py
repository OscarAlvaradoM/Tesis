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
    
    
        Métodos:
        constructor(nvx,delta): inicia los atributos nvx y delta
        destructor(): delete atributes
        alloc(n): asigna arreglos con ceros a los atributos de coeficientes (aP, aE, etc.) con
                  el objetivo de reservar memoria
        setVolumes(nvx): set atributo nvx
        setDelta(delta): set atributo delta
        aP():get aP
        aW(): get aW
        aWW(): get aWW
        aE(): get aE
        aEE(): get aEE
        sU(): get sU                
        bcDirichlet(wall,phi): ajusta los coeficientes de la frontera 'wall' (puede ser 'LEFT_WALL' o 'RIGHT_WALL')
                               de acuerdo a la condición de frontera con valor 'phi'
        bcNeumman(wall,flux): ajusta los coeficientes de la frontera 'wall' (puede ser 'LEFT_WALL' o 'RIGHT_WALL')
                               de acuerdo a la condición de frontera con flujo de valor 'flux'
        setsU(q): set atributo Su
        setSp(Sp): set atributo Sp

        
    Attributes:
        dim: number of dimensions in the physical domain (int)
        mesh: mesh defined to discretize the physical domain ()
        N: total number of nodes (int)
        nvzyx: tuple with the number of volumes in z, y and x direction in that order (tuple of int) 
        Su: array of the independet coefficients, asociated with sources, for the discrete equations of the volumes (numpy array)
        Sp: array of the dependent coefficients, asociated with sources, for the discrete equations of the volumes (numpy array)
        aP: array with all the coefficients corresponding to the central node in the discrete equations (numpy array)
        aE: array with all the coefficients corresponding to the east node in the discrete equations (numpy array)
        aW: array with all the coefficients corresponding to the west node in the discrete equations (numpy array)
        aN: array with all the coefficients corresponding to the north node in the discrete equations (numpy array)
        aS: array with all the coefficients corresponding to the south node in the discrete equations (numpy array)
        aT: array with all the coefficients corresponding to the top node in the discrete equations (numpy array)
        aB: array with all the coefficients corresponding to the bottom node in the discrete equations (numpy array)
        
    """    

    def __init__(self, mesh):
        
        self._dim = mesh.dim()
        self._mesh = mesh
        self._N = ( mesh.nvx() * mesh.nvy() * mesh.nvz() )
        self._nvzyx = mesh.nvzyx()
        self._Su = np.zeros(self._nvzyx)  ;  self._Sp = None
        self._aN = None ; self._aS = None ; self._aT = None ; self._aB = None
        self.alloc_a()
        
    def alloc_a(self):
        
        """ Depending on the dimension of the physical domain, this method initialize
        the coefficients as 3D arrays full of zeros"""
        
        nvzyx = self._nvzyx
        self._aP = np.zeros(nvzyx)
        dim = self._dim
        self._aE = np.zeros(nvzyx)
        self._aW = np.zeros(nvzyx)
        if dim > 1:
            self._aN = np.zeros(nvzyx)
            self._aS = np.zeros(nvzyx)
        if dim == 3:
            self._aT = np.zeros(nvzyx)
            self._aB = np.zeros(nvzyx)
              

    def setSu(self, su):
        """
        Set the 'Su' attribute as a 3D array by taking the parameter 'su'
        given and multiplying it with the volumes volumes. The parameter 'su'
        can be a float or a numpy array. If an array is given it should have the
        same shape that the 'mesh.vols()', that is (nvz,nvy,nvx), and if a float
        is given then the 3D array 'Su' attribute is stored having a constant
        value for all it's elements.
        
        su: independet coefficients, asociated with sources, for the discrete equations of the nodes (float or numpy array)
        """        
        mesh = self._mesh
        vols = mesh.vols()
        self._Su = vols * su
        
    def saveAsSu(self, su):
        """

        """        
        self._Su = su
        
    def saveSu(self, su):
        """

        """        
        self._Su = su

    def setSp(self, sp):
        """
        Set the 'Sp' attribute as a 3D array by taking the parameter 'sp'
        given and multiplying it with the volumes volumes. The parameter 'sp'
        can be a float or a numpy array. If an array is given it should have the
        same shape that the 'mesh.vols()', that is (nvz,nvy,nvx), and if a float
        is given then the 3D array 'Sp' attribute is stored having a constant
        value for all it's elements.
        
        sp: dependent coefficients, asociated with sources, for the discrete equations of the nodes (float or numpy array)
        """
        
        mesh = self._mesh
        ####-----------No es necesario guardar el atributo Sp ----#######
        ## Para ahorrar memoria se puede corregir aP sin guardar Sp en memoria ##
        vols = mesh.vols()
        self._Sp = vols * sp
        self._aP -= self._Sp

        
    def setaP(self, aP):
        self._aP = aP
    
    def setaE(self, aE):
        self._aE = aE

    def setaW(self, aW):
        self._aW = aW
        #self._aE = aW

    def setaN(self, aN):
        self._aN = aN

    def setaS(self, aS):
        self._aS = aS

    def setaT(self, aT):
        self._aT = aT

    def setaB(self, aB):
        self._aB = aB
        
    def setDiffusion(self, Gamma):
        """This method makes an instance of the Diffusion class and uses it to
        update the aP,aE,aW, ... atributes.
        
        Gamma: The diffusion coeffcient of the transport equation, if the diffusion
               coefficient is constant then a float can be given but in other cases
               a function of the coordinates, denoted as Gamma=f(x,y,z) should be 
               used. (float or function)"""
        
        dim = self._dim
        malla = self._mesh
        D = Diffusion(malla, Gamma)
        DE = D.DE()
        DW = D.DW()        
        self._aE -= DE
        self._aW -= DW
        self._aP +=  DE + DW
        if dim > 1:
            DN = D.DN()
            DS = D.DS()
            self._aN -= DN
            self._aS -= DS
            self._aP += DN + DS            
        if dim == 3:
            DT = D.DT()
            DB = D.DB()
            self._aT -= DT
            self._aB -= DB
            self._aP += DT + DB   
            
            
    def under(self, alpha, velComp):
        self._aP = self._aP / alpha
        self._Su += (1. - alpha) * self._aP * velComp
        
            
    def setTime(self, phi_old, dt, rho,scheme = 'implicit'):
        
        if scheme == 'implicit':
        
            mesh = self._mesh
            vols = mesh.vols()
            rv_t = (rho * vols) / dt
            self._aP += rv_t
            self._Su += phi_old * rv_t  ### Este paso podría sacarse de este método y meterse a un nuevo método para separar la modificacion a aP y a Su y así no tener que redefinir los coeficientes en cada paso temporal

    def setTime2(self, dt, rho, scheme = 'implicit'):
        
        if scheme == 'implicit':
        
            mesh = self._mesh
            vols = mesh.vols()
            rv_t = (rho * vols) / dt
            self._aP += rv_t
            
    def doStepTime(self, phi_old, dt, rho, scheme = 'implicit'):
        
        if scheme == 'implicit':
        
            mesh = self._mesh
            vols = mesh.vols()
            rv_t = (rho * vols) / dt
            self._Su += phi_old * rv_t 
            
    def Su(self):
        return self._Su
    
    def Sp(self):
        return self._Sp
    
    def N(self):
        return self._N
    
    def aP(self):
        return self._aP
    
    def aE(self):
        return self._aE
    
    def aW(self):
        return self._aW

    def aN(self):
        return self._aN

    def aS(self):
        return self._aS

    def aT(self):
        return self._aT

    def aB(self):
        return self._aB
    
    def mesh(self):
        return self._mesh
    
    def nvzyx(self):
        return self._nvzyx
    
    