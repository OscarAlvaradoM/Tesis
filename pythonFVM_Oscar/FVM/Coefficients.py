import numpy as np
from Diffusion import Diffusion
from Convection import Convection
class Coefficients():
    """
    Class that defines the coefficients that every node needs to form it's
    discrete equation that represent the diferential equation. All the coefficients
    needed for a scalar quantity (also called variable or unkown) are contained
    in one object (instance of this class).
    """    

    def __init__(self, mesh):
        self.dim = mesh.dim
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
        self.aP = np.zeros(vols)
        dim = self.dim
        self.aE = np.zeros(vols)
        self.aW = np.zeros(vols)
        self.Sp = np.zeros(vols)
        self.Su = np.zeros(vols)

        self.aN = np.zeros(vols)
        self.aS = np.zeros(vols)

        self.aT = np.zeros(vols)
        self.aB = np.zeros(vols)

    def set_diffusion(self, gamma):
        """This method makes an instance of the Diffusion class and uses it to
        update the aP,aE,aW, ... atributes.
        
        Gamma: The diffusion coeffcient of the transport equation, if the diffusion
               coefficient is constant then a float can be given but in other cases
               a function of the coordinates, denoted as Gamma=f(x,y,z) should be 
               used. (float or function)"""
        
        dim = self.dim
        malla = self._mesh
        diffusion = Diffusion(malla, gamma)
        east_diff, sp_e, su_e = diffusion.east_diffusion()
        west_diff, sp_w, su_w = diffusion.west_diffusion()
        self.aE -= east_diff
        self.aW -= west_diff
        self.Sp -= sp_e + sp_w
        self.Su += su_e + su_w
        self.aP -=  self.aE + self.aW
        if dim > 1:
            north_diff, sp_n, su_n = diffusion.north_diffusion()
            south_diff, sp_s, su_s = diffusion.south_diffusion()
            self.aN -= north_diff
            self.aS -= south_diff
            self.Sp -= sp_n + sp_s
            self.Su += su_n + su_s
            self.aP -= self.aN + self.aS
        if dim == 3:
            top_diff, sp_t, su_t = diffusion.top_diffussion()
            bottom_diff, sp_b, su_b = diffusion.bottom_diffusion()
            self.aT -= top_diff
            self.aB -= bottom_diff
            self.Sp -= sp_t + sp_b
            self.Su += su_t + su_b
            self.aP -= top_diff + bottom_diff
        self.aP -= self.Sp
        
    def set_convection(self, rho, vel, scheme='central'):
        """
        This method makes an instance of the convection class and uses it to
        update the aP, aE, aW, ... atributes.
        
        rho: float or 3D array that represents the density, if a float is used
             it is assumed that the density is constant in all the physical 
             domain. (float/numpy array)
        
        vel: tuple of 3D arrays that represents the velocity, each 3D array in
             the tuple represent a component of the velocity. It can be denoted 
             as vel=(u,v,w) where u, v and w are 3D arrays with values for all
             of the volume borders positions. (tuple of numpy arrays) 
        
        scheme: scheme to be used for the convection terms. Strings 'central' and 
                'upwind1' are currently valid (string).
        """
        
        dim = self.dim
        malla = self._mesh
        convection = Convection(malla, rho, vel, scheme)
        east_conv, sp_e, su_e = convection.east_convection()
        west_conv, sp_w, su_w = convection.west_convection()
        east_f = convection.get_f(east_conv, "E")
        west_f = convection.get_f(west_conv, "W")
        self.aE += east_f
        self.aW -= west_f
        self.aP -= east_f + sp_e - sp_w - west_f
        self.Su += su_e + su_w
        self.Sp -= -sp_e + sp_w
        if dim > 1:
            north_conv, sp_n, su_n = convection.south_convection()
            south_conv, sp_s, su_s = convection.south_convection()
            north_f = convection.get_f(north_conv, "N")
            south_f = convection.get_f(south_conv, "S")
            self.aE += north_f
            self.aW -= south_f
            self.aP -= north_f + sp_n - sp_s - south_f
            self.Su += su_n + su_s
            self.Sp -= -sp_n + sp_s
        if dim == 3:
            top_conv, sp_t, su_t = convection.top_convection()
            bottom_conv, sp_b, su_b = convection.bottom_convection()
            top_f = convection.get_f(top_conv, "T")
            bottom_f = convection.get_f(bottom_conv, "B")
            self.aE += top_f
            self.aW -= bottom_f
            self.aP -= top_f + sp_t - sp_b - bottom_f
            self.Su += su_t + su_b
            self.Sp -= -sp_t + sp_b
            
    def get_Su(self):
        return self.Su
    
    def get_Sp(self):
        return self.Sp
    
    def get_aP(self):
        return self.aP
    
    def get_aE(self):
        return self.aE
    
    def get_aW(self):
        return self.aW

    def get_aN(self):
        return self.aN

    def get_aS(self):
        return self.aS

    def get_aT(self):
        return self.aT

    def get_aB(self):
        return self.aB
