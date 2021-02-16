#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
@author: José ft. Óscar
"""

import numpy as np
import plotly.graph_objects as go

class Mesh():
    """
    Class that defines the scalar mesh in which FVM will be apllied. This class contains methods for defining
    borders, getting volumes/nodes positions, visualization of mesh, etc. The term volumes refers to points, 
    more specifically, the term volumes is used for the center of the volumes but the border nodes are not
    included. The terms nodes acounts fot all the points (interior and border) of the mesh. For every node
    there is a tag; the tags is an string beginning with I,D,N or S indicating whether the node is interior,
    dirichlet, neumann or source, accordingly. The atributes that are lists or arrays that contains an element 
    for every node, like the 'tags' atribute, are sorted sweeping in the X, then Y an finally Z diretion;
    e.g. tags=[tag1,tag2,tag3,tag4,tag5,...] corresponds to nodes coordinates sorted like [(1,1,1),(2,1,1),(3,1,1),(1,2,1),(2,2,1),...]
    
    
    Methods:
        constructor(dim,volumes,lenght): instance a mesh
        totalDomNodes(): get the total number of nodes in the domain
        uniformGrid(): set mesh atributes assuming the mesh is uniform (it might be used by the constructor)
        uniformVector(li,nvi): returns the 'nvi' volumes postions, nodes positions and nodes separation for a line segment of lengh 'li'
        setDominio(domino): set mesh atributes according to the nodes positions given by parameter 'dominio'
        setDetlas(dominio): Obtain the separation between the positions if the array 'dominio'
        totalDomNodes(): Get the total number of nodes
        tagsWallsOff(): Define a tag Off for all the nodes not used in the defined dimension
        NdomX(): get the number of nodes in the X direction
        NdomY(): get the number of nodes in the Y direction
        NdomZ(): get the number of nodes in the Z direction
        info(): prints some relevant information about the mesh
        createCoord(): create and stores as atributtes the nodes coordinates
        delCoord(): delete the coordX,coordY and coordZ attributes
        getCoord(): get the nodes coordinates
        coordList(): get the nodes coordinates in a printable fashion
        coordTags(): get the nodes coordinates and their respective tags
        printCoordTags(includeOffs): get the nodes coordinates and their respective tags and prints them in a convenient way
        
    Attributes:
        nvx: number of volumes in dimension X (int)
        nvy: number of volumes in dimension Y (int)
        nvz: number of volumes in dimension Z (int)
        lx: length of the domain in the X direction (float)
        ly: length of the domain in the Y direction (float)
        lz: length of the domain in the Z direction (float)
        X: list of 'x' coordinates, without repetition, of volumes (list of floats or numpy array)
        Y: list of 'y' coordinates, without repetition, of volumes (list of floats or numpy array)
        Z: list of 'z' coordinates, without repetition, of volumes (list of floats or numpy array)
        dominioX: list of 'x' coordinates, without repetition, of nodes (list of floats or numpy array)
        dominioY: list of 'y' coordinates, without repetition, of nodes (list of floats or numpy array)
        dominioZ: list of 'z' coordinates, without repetition, of nodes (list of floats or numpy array)
        deltaX: list of separations between the 'dominioX' elements (numpy array)
        deltaY: list of separations between the 'dominioY' elements (numpy array)
        deltaZ: list of separations between the 'dominioZ' elements (numpy array)
        coordX: sorted list of 'x' coordinates, with repetition, of nodes (numpy array)
        coordY: sorted list of 'y' coordinates, with repetition, of nodes (numpy array)
        coordZ: sorted list of 'z' coordinates, with repetition, of nodes (numpy array)
        dirichTagDict: dictionary representing the relation between the dirichlet tags and their values (dictionary)
        neumTagDict: dictionary representing the relation between the neumann tags and their values (dictionary)
        sourcTagDict: dictionary representing the relation between the source tags and their values (dictionary)
        dim: number of dimensions in the physical domain (int)
        volumes: number of volumes in each dimension (tuple of ints or int)
        lengths: lenght of domain in each dimension (tuple of ints or int)
        tags: sorted list that has all the nodes tags (list of strings)
        autoMesh: flag to indicate whether coordinates were calculated automatically using the volumes and lengths attributes (True/False)
        emptyCoord: flag to indicate whether the coordinates are left unstored or not (True/False)
        pressed: flag that becomes True when the mesh, and its properties, wont require any changes (True/False)
        intWallVols:
        hasWallW: list of all the nodes index; 0,1,2,3,...,etc, that have a Wall in their West node (list)
        hasWallE: list of all the nodes index, that have a Wall in their East node (list)
        hasWallN: list of all the nodes index, that have a Wall in their North node, the list is empty if dim=1 (list)
        hasWallS: list of all the nodes index, that have a Wall in their South node, the list is empty if dim=1 (list)
        hasWallT: list of all the nodes index, that have a Wall in their Top node, the list is empty if dim<3 (list)
        hasWallB: list of all the nodes index, that have a Wall in their Bottom node, the list is empty if dim<3 (list)
        
    """
    
    def __init__(self, dim, volumes = None, lengths = None):
        """
        Constructor of the Mesh class. If the parameters 'volumes' and 'lengths' are given then 
        the instance atributes are defined assuming the mesh is uniform via the uniformGrid() method.
        When volumes or lenghts are leaved as None then the domain nodes must be defined 
        using setDmonio() just after making the instance of the class.

        dim: number of dimensions in the physical domain (int)
        volumes: number of volumes in each dimension (tuple of ints or int)
        lengths: lenght of domain in each dimension (tuple of ints or int)
        """

        self.__volumes = (1,1,1)
        self.__lengths = (0.01, 0.01, 0.01)

        #----default values for positions and separations of the grid-------------
        self.__coords = [(0.005,) for _ in range(3)] # Coordenadas de los centros de los volúmenes
        self.__dominios = [tuple([0.]+[self.__coords[i][0]]+[self.__lengths[i]]) for i in range(3)]
        self.__deltas = [(0,) for _ in range(3)]
        # Hasta aquí tenemos un cubito

        self.__tags = {} # El etiquetado de todos los nodos sin fronteras
        self.__tags_fronteras = {} # El etiquetado de las fornteras

        self.__autoMesh = True      # Al parecer es una bandera para calcular la posición de los nodos atmte
        self.__emptyCoord = True    # Nos dice si las coordenadas del dominio (self.__coords_dom) están guardadas
        self.__pressed = False      # No sé

        self.__intWallNodes = None  # No sé
        
        self.__dim = dim
        #---if a parameter is not a tuple, that parameter is transformed into a tuple
        if isinstance(volumes, int):  self.__volumes = (volumes, 1, 1)
        if isinstance(lengths, (int, float)):  self.__lengths = (lengths, lengths/10, lengths/10)
        
        # Si los parámetros son tuplas (pero no necesariamente sería una tupla de 3, arreglamos eso:
        if isinstance(volumes, tuple):
            faltan = 3 - len(volumes)
            self.__volumes = volumes + tuple([1 for i in range(faltan)])
        if isinstance(lengths, tuple): 
            faltan = 3 - len(lengths)
            self.__lengths = lengths + tuple([lengths[0]/10 for i in range(faltan)])
        #---------------------------------------------------------------
        
        # if volumes and lengths are given, initialize values acording to dimension
        if (volumes and lengths):       
            self.uniformGrid()
            self.init_tags()
            self.init_tags_fronteras()
    
    
    def uniformGrid(self):
        l = np.array(self.__lengths) # Para el manejo con numpy
        v = np.array(self.__volumes)
        d = l/v # Separación entre todos los nodos de cada dimensión
        start = d/2 # La frontera "inicial" del arreglo
        stop = l-d/2 # La frontera "final" del arreglo
        self.__coords = [tuple(np.linspace(strt, stp, vol)) for strt, stp, vol in list(zip(start, stop, v))] # Meshgrid posible
        dominios = [np.insert(arr,(0,len(arr)),[0, l[idx]]) for idx, arr in enumerate(self.__coords)] # Coordenadas + fronteras
        self.__deltas = [self.setDeltas(dom)  if len(self.setDeltas(dom)) != 0 else (dom[-1],) for dom in dominios]  # Separación entre los nodos (Aquí hay que ver cómo es cuando tenemos un grid de 2x1x1 ya cuando se haga el FVM
        self.__dominios = [tuple(dom) for dom in dominios]
        self.__faces = [tuple((np.array(coords[:-1]) + np.array(coords[1:]))/2) for coords in self.__coords]
    
    
    def setDeltas(self, dominio):
        return tuple((dominio[1:]-dominio[:-1])[1:-1])
    
    
    # Creo que esto no se usará, pero estuvo chida la deducción, lo dejo de todos modos xd
    #def totalDomNodes(self):
    #    d_1 = 6*self.__volumes[0] + 1
    #    d_2 = self.__volumes[1]*d_1 - self.__volumes[0]*(self.__volumes[1] - 1)
    #    d_3 = self.__volumes[2]*d_2 - self.__volumes[0]*self.__volumes[1]*(self.__volumes[2] - 1)
    #    return d_3
    
    def init_tags_fronteras(self):
        """
        Método para etiquetar las fronteras dependiendo de la dimensión, sólo se les da la propiedad de existir o no
        existir.
        """
        X, Y, Z = [len(dom) for dom in self.__dominios]
        for z in range(Z):
            for y in range(Y):
                for x in range(X):
                    t = b = n = s = "Off"
                    e = w = "ON"
                    if self.__dim > 1: 
                        n = s = "ON"
                        if self.__dim == 3: t = b = "ON"
                    # El siguiente cacho de código es para saber si nos encontramos con una frontera
                    if x==0 or y==0 or z==0 or x==(X-1) or y==(Y-1) or z==(Z-1):
                        var = None
                        if y != 0 and y != (Y - 1):
                            if z != 0 and z != (Z - 1):
                                if x == 0: var = "W"; value = w
                                elif x == (X - 1): var = "E"; value = e
                                else: continue
                            elif x != 0 and x != (X - 1):
                                if z == 0: var = "B"; value = b
                                elif z == (Z - 1): var = "T"; value = t
                                else: continue
                            else: continue
                            self.__tags_fronteras[f"{x}{y}{z}"] = {"frontera": {var: value},
                                                 "coord": [self.__dominios[0][x], self.__dominios[1][y], self.__dominios[2][z]],
                                                                  "cond": {}} 
                        elif z != 0 and z != (Z - 1):
                            if x != 0 and x != (X - 1):
                                if y == 0: var = "S"; value = s
                                elif y == (Y - 1) : var = "N"; value = n
                                self.__tags_fronteras[f"{x}{y}{z}"] = {"frontera": {var: value},
                                                 "coord": [self.__dominios[0][x], self.__dominios[1][y], self.__dominios[2][z]],
                                                                      "cond": {}} 
                        else: continue
    
    def init_tags(self):
        """
        Método que etiqueta las caras adyacentes de cada volumen dependiendo de la geometría. Pone un 0 (cero) cuando es una 
        frontera, una 'F' cuando es una cara interna y un 'Off' cuando no se está contando esa cara por las dimensiones del 
        problema. 
        """
        X, Y, Z = self.__volumes
        for z in range(1,Z+1):
            for y in range(1,Y+1):
                for x in range(1,X+1):                   
                    t = b = n = s = "Off"
                    e = w = "F"
                    if x == 1: w = {}
                    elif x == X: e = {}
                    if self.__dim > 1:
                        n = s = "F"
                        if y == 1: s = {}
                        elif y == Y: n = {}
                        if self.__dim == 3:
                            t = b = "F"
                            if z == 1: b = {}
                            elif z == Z: t = {}

                    self.__tags[f"{x}{y}{z}"] = {"E": e, "W": w, "N": n, "S": s, "T": t, "B": b, 
                                             "coord": [self.__dominios[0][x], self.__dominios[1][y], self.__dominios[2][z]]}
            
    def tag_wall(self, direction, tag, value):
        """
        Método para etiquetar fronteras dada la dirección, el tipo de condición de frontera y el valor.
        """
        for key in self.__tags.keys():
            if isinstance(self.__tags[key][direction], dict):
                self.__tags[key][direction][tag] = value
        for key in self.__tags_fronteras.keys():
            if self.__tags_fronteras[key]["frontera"].get(direction) == "ON":
                self.__tags_fronteras[key]["cond"][tag] = value

    def tag_wall_dirichlet(self, direction, value, coords=None):
        """
        Método para etiquetar fronteras con condición de Dirichlet dados ciertos valores.
        """
        if coords:
            for idx, key in enumerate(coords):
                if key in list(self.__tags.keys()):
                    self.__tags[key][direction[idx]]["D"] = value[idx]
                elif key in list(self.__tags_fronteras.keys()):
                    self.__tags_fronteras[key]["cond"]["D"] =  value[idx]
        else:
            if isinstance(direction, list):
                for idx, direct in enumerate(direction):
                    self.tag_wall(direct, "D", value[idx])
            else:
                self.tag_wall(direction, "D", value)
                    
    def tag_wall_neumann(self, direction, value, coords=None):
        """
        Método para etiquetar fronteras con condición de Neumann dados ciertos valores.
        """
        if coords:
            for idx, key in enumerate(coords):
                if key in list(self.__tags.keys()):
                    self.__tags[key][direction[idx]]["N"] = value[idx]
                elif key in list(self.__tags_fronteras.keys()):
                    self.__tags_fronteras[key]["cond"]["N"] = value[idx]
        else:
            if isinstance(direction, list):
                for idx, direct in enumerate(direction):
                    self.tag_wall(direct, "N", value[idx])
            else:
                self.tag_wall(direction, "N", value)
                    
    def tag_wall_source(self, direction, value, coords=None):
        """
        Método para etiquetar fronteras con condición de Neumann dados ciertos valores.
        """
        if coords:
            for idx, key in enumerate(coords):
                if key in list(self.__tags.keys()):
                    self.__tags[key][direction[idx]]["S"] = value[idx]
                elif key in list(self.__tags_fronteras.keys()):
                    self.__tags_fronteras[key]["cond"]["S"] = value[idx]
        else:
            if isinstance(direction, list):
                for idx, direct in enumerate(direction):
                    self.tag_wall(direct, "S", value[idx])
            else:
                self.tag_wall(direction, "S", value)
                
    
    def set_dominio(self, dominio, faces=None):
        """
        Método para definir el dominio de estudio dadas unas coordenadas en forma de tupla.
        """
        
        self.__autoMesh = False # La posición de los nodos no se calcula en automático
        
        # Si 'dominio' no es tupla, transforma 'dominio' a la tupla unidimensional (dominio,)
        # Tendría que ser una tupla de tuplas/listas/arreglos para que sea válido.
        if not isinstance(dominio, (tuple, int, float)): # Creo que si es una lista o un arreglo, no funciona enteros o float
            tupla = (tuple(dominio), self.__dominios[1], self.__dominios[2])
            dominio = tupla
        # Asigna los atributos de la mesh correspondientes    
        self.__dominios = [tuple(dominio[i]) for i in range(3)]
        self.__coords = [tuple(dominio[i][1:-1]) for i in range(3)]
        self.__lengths = tuple([dominio[i][-1] for i in range(3)])
        self.__volumes = tuple([len(dominio[i][1:-1]) for i in range(3)])
        self.__deltas = [self.setDeltas(np.array(dominio[i])) for i in range(3)]
        
        if faces: self.__faces = faces
        else: self.__faces = [tuple((np.array(coords[:-1]) + np.array(coords[1:]))/2) for coords in self.__coords]
            
        self.init_tags()
        self.init_tags_fronteras()
    
    
    def info(self):
        """
        Método para imprimir información relevante del mallado
        """
        print('=====================================')
        print('     MESH INFORMATION   ')
        print('=====================================')
        print("\nMesh type: Cartesian")
        print(f"Number of dimensions of mesh: {self.__dim}")
        variables = "X Y Z".split()
        for idx in range(self.__dim):
            var = variables[idx]
            print(f"\n ----- {var}-axis -----")
            print(f"Number of {var.lower()} volumes: {self.__volumes[idx]}")
            print(f"Lenght {var.lower()} of problem domain: {self.__lengths[idx]}")
            print(f"List of {var.lower()} positions of volumes: \n{self.__coords[idx]}")
            print(f"List of {var.lower()} positions of domain nodes: \n{self.__dominios[idx]}")
            
            
    def draw(self):
        """
        Método para graficar la malla. Este método se invoca hasta que se hayan inizializado todas las condiciones
        de frontera.
        """
        # Graficamos las fronteras, sean o no activas
        dic_colors = {"D": "gold", "N": "red", "S": "magenta", "Off": "gray"}
        condiciones = [list(self.__tags_fronteras[key]["cond"].keys())[0] if list(self.__tags_fronteras[key]["frontera"].values())[0] == "ON" else "Off" for key in list(self.__tags_fronteras.keys())]
        colores = [dic_colors[cond] for cond in condiciones]
        # Obtenemos las coordenadas de las fronteras y de los nodos internos.
        coordenadas = [] # Aquí se pondrán las coordenadas de las fornteras
        coord = [] # Aquí las coordendas de los nodos internos
        for i in range(3):
            coordenadas.append([self.__tags_fronteras[key]["coord"][i] for key in list(self.__tags_fronteras.keys())])
            coord.append([self.__tags[key]["coord"][i] for key in list(self.__tags.keys())])
        fig = go.Figure(data = go.Scatter3d(x = coordenadas[0], y = coordenadas[1], z = coordenadas[2],
                                              mode = 'markers', marker = dict(color = colores, symbol = "square", size = 2)))
        fig.add_trace(go.Scatter3d(x = coord[0], y = coord[1], z = coord[2],
                                              mode = 'markers', marker = dict(color = "blue", size = 5)))
        fig.show()
    
    
    def u_stag_def(self):
        """
        Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the X axis.
        """
        stgX = np.zeros(self.__volumes[0]+1)
        stgX[-1] = self.__lengths[0]
        stgX[1:-1] = 0.5*(self.X()[:-1] + self.X()[1:])
        return (stgX,self.dominioY(),self.dominioZ())
    
    
    def v_stag_def(self):
        """
        Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the Y axis.
        """
        
        stgY=np.zeros(self.__volumes[1]+1)
        stgY[-1]=self.__lengths[1]
        stgY[1:-1]=0.5*(self.Y()[:-1]+self.Y()[1:])
        return (self.dominioX(),stgY,self.dominioZ())
    
    
    def w_stag_def(self):
        """
        Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the Z axis.
        """
        
        stgZ=np.zeros(self.__volumes[2]+1)
        stgZ[-1]=self.__lengths[2]
        stgZ[1:-1]=0.5*(self.Z()[:-1]+self.Z()[1:])
        return (self.dominioX(),self.dominioY(),stgZ)


    def vols(self):
        """
        Returns a 3-dimensional numpy array containing all the value of the volumes,
         that is lenght times width times height, for each one.
         """
        dim = self.__dim
        dx=self.deltaX() ; dy=self.deltaY() ; dz=self.deltaZ()
        dx_vol = (dx[1:]+dx[:-1])/2. 
        dy_vol = (dy[1:]+dy[:-1])/2.
        dz_vol = (dz[1:]+dz[:-1])/2.
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        vols_dim_grid = np.meshgrid(dy_vol,dz_vol,dx_vol)
        return vols_dim_grid[0]*vols_dim_grid[1]*vols_dim_grid[2]
    
    def areas_x(self):
        """
        Returns a 3-dimensional numpy array containing the areas of the faces in
        the X direction. The numpy array contains an area for each volume. 
        """
        
        dim = self.__dim
        dy=self.deltaY() ; dz=self.deltaZ()
        nvx=self.nvx()
        dy_vol = (dy[1:]+dy[:-1])/2.  #longitud y (ancho) del volumen de control
        dz_vol = (dz[1:]+dz[:-1])/2.  #longitud z (alto) del volumen de control
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        areas_grid =np.meshgrid(dy_vol,dz_vol,np.ones(nvx) )
        return areas_grid[0]*areas_grid[1]

    def areas_y(self):
        """
        Returns a 3-dimensional numpy array containing the areas of the faces in
        the Y direction. The numpy array contains an area for each volume. 
        """

        dim = self.__dim
        dx=self.deltaX() ; dz=self.deltaZ()
        nvy=self.nvy()
        dx_vol = (dx[1:]+dx[:-1])/2.  #longitud y (ancho) del volumen de control
        dz_vol = (dz[1:]+dz[:-1])/2.  #longitud z (alto) del volumen de control
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        areas_grid =np.meshgrid(np.ones(nvy),dz_vol, dx_vol)
        return areas_grid[2]*areas_grid[1]

    def areas_z(self):
        """
        Returns a 3-dimensional numpy array containing the areas of the faces in
        the Z direction. The numpy array contains an area for each volume. 
        """  
        
        dim = self.__dim
        dx=self.deltaX() ; dy=self.deltaY()
        nvz=self.nvz()
        dx_vol = (dx[1:]+dx[:-1])/2.  #longitud y (ancho) del volumen de control
        dy_vol = (dy[1:]+dy[:-1])/2.  #longitud y (ancho) del volumen de control
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        areas_grid =np.meshgrid(dy_vol,np.ones(nvz),dx_vol )
        return areas_grid[2]*areas_grid[0]