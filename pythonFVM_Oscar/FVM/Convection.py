import numpy as np
class Convection():
    """
    rho: float or 3D array that represents the density, if a float is used
         it is assumed that the density is constant in all the physical 
         domain. (float/numpy array)

    vel: tuple of 3D arrays that represents the velocity, each 3D array in
         the tuple represent a component of the velocity. It can be denoted 
         as vel=(u,v,w) where u, v and w are 3D arrays with values for all
         of the volume borders positions. (tuple of numpy arrays) 

    scheme: scheme to be used for the advection terms. Strings 'central' and 
            'upwind1' are currently valid (string).
    """
    def __init__(self,mesh,ρ,vel,scheme):
        self._mesh   = mesh
        self._dim    = mesh.dim
        self._vel    = vel
        self._scheme = scheme
        
        self._ρ_constant = 0.
        self._ρ = None
        
        if isinstance(ρ, (int,float)):
            self._ρ_constant = ρ
            self._ρ = self.funcion_const
        # Se asume que es una función en este caso
        else:
            self._ρ = ρ
            
    def funcion_const(self, x, y, z):
        return self._ρ_constant
        
    def get_convection_coef(self, direction):
        idx, idx_1 = 0,0
        if direction in ["E", "N", "T"]: idx_1=-1
        if direction in ["N", "S"]: idx=1
        elif direction in ["T", "B"]: idx = 2
        mesh = self._mesh
        faces = mesh.faces[idx]
        x, y, z = np.meshgrid(mesh.coords[0], mesh.coords[1], mesh.coords[2])

        # Para obtener los coeficientes
        λ_1 = lambda x, d: x[1:] if d in ["E","N","T"] else x[:-1]
        coord = [var if i!=idx else λ_1(faces, direction) for (i, var) in enumerate([x,y,z])]
        ρ = self._ρ(coord[0], coord[1], coord[2])
        areas = mesh.get_area(direction=idx)
        source = areas.copy()
        coord_2 = [slice(None,None) if var!=idx else idx_1 for var in range(3)]
        areas[coord_2[0],coord_2[1],coord_2[2]] = 0
        λ_2 = lambda d: slice(None,-1) if d in ["E","N","T"] else slice(1,None)
        coord_3 = [slice(None,None) if var!=idx else λ_2(direction) for var in range(3)]
        source[coord_3[0],coord_3[1],coord_3[2]] = 0
        velocity = self._vel[idx][coord_3[0],coord_3[1],coord_3[2]]
        conv = ρ*velocity*areas
        
        # Aquí obtenemos Sp
        sp = ρ*velocity*source
        condicion = self.get_mask_boundaries_Sp(mesh, direction)
        λ_3 = lambda d: slice(-1,None) if d in ["E","N","T"] else slice(None,1)
        coord_4 = [slice(None,None) if var!=idx else λ_3(direction) for var in range(3)]
        tmp_sp = sp[coord_4[0],coord_4[1],coord_4[2]]
        sp[coord_4[0],coord_4[1],coord_4[2]] = tmp_sp*np.array(condicion).reshape(tmp_sp.shape)
        
        # Aquí obtenemos Su
        conds_su = self.get_mask_boundaries_Su(mesh, direction, gamma=self._ρ)
        su = ρ*source*velocity
        tmp_su = su[coord_4[0],coord_4[1],coord_4[2]]
        su[coord_4[0],coord_4[1],coord_4[2]] = tmp_su*np.array(conds_su).reshape(tmp_su.shape)
        return conv, sp, su

    def east_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aE coefficient 
        caused by the convection effect.
        """
        east_d, sp, su = self.get_convection_coef("E")
        return east_d, sp, su

    def west_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aW coefficient
        caused by the convection effect.
        """        
        west_d, sp, su = self.get_convection_coef("W")
        return west_d, sp, su

    def north_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aN coefficient
        caused by the convection effect.
        """
        north_d, sp, su = self.get_convection_coef("N")
        return north_d, sp, su

    def south_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aS coefficient
        caused by the convection effect.
        """
        south_d, sp, su = self.get_convection_coef("S")
        return south_d, sp, su

    def top_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aT coefficient
        caused by the convection effect.
        """
        top_d, sp, su = self.get_convection_coef("T")
        return top_d, sp, su

    def bottom_convection(self):
        """
        Gives a 3D numpy array that should be used to correct the aB coefficient
        caused by the convection effect.
        """
        bottom_d, sp, su = self.get_convection_coef("B")
        return bottom_d, sp, su
    
    
    def get_f(self, coef, sentido="sup"):
        sign = 1
        if sentido == "inf": sign = -1
        
        scheme = self._scheme
        if scheme == 'central':
            f = sign*0.5*coef
        if scheme == 'upwind1':
            ceros = np.zeros(coef.shape)
            f = -np.maximum(-sign*coef, ceros)
        return f
    
    def get_mask_boundaries_Sp(self, malla, direction):
        tags_fronteras = malla._Mesh__tags_fronteras
        condicion = []
        dict_cond = {"I":0, "N":0, "D":1}
        for tag in tags_fronteras:
            if list(tags_fronteras[tag]["frontera"].keys())[0] == direction:
                cond = list(tags_fronteras[tag]["cond"].keys())[0]
                condicion.append(dict_cond[cond])
        return condicion
    
    def get_mask_boundaries_Su(self, malla, direction, gamma):
        tags_fronteras = malla._Mesh__tags_fronteras
        condicion = []
        dict_cond = {"I":0, "N":0, "D":1}
        for tag in tags_fronteras:
            if list(tags_fronteras[tag]["frontera"].keys())[0] == direction:
                cond = list(tags_fronteras[tag]["cond"].keys())[0]
                if cond == "I":
                    condicion.append(0)
                elif cond == "N":
                    condicion.append(tags_fronteras[tag]["cond"][cond])
                elif cond == "D":
                    x,y,z = tags_fronteras[tag]["coord"]
                    condicion.append(tags_fronteras[tag]["cond"][cond]*gamma(x,y,z))
        return condicion
