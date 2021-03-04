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
    """
    
    def __init__(self, dim, volumes = None, lengths = None):
        """
        Constructor of the Mesh class. If the parameters 'volumes' and 'lengths' are given then 
        the instance atributes are defined assuming the mesh is uniform via the uniform_grid() method.
        When volumes or lenghts are leaved as None then the domain nodes must be defined 
        using setDmonio() just after making the instance of the class.

        dim: number of dimensions in the physical domain (int)
        volumes: number of volumes in each dimension (tuple of ints or int)
        lengths: lenght of domain in each dimension (tuple of ints or int)
        """

        self.volumes = (1,1,1)
        self.lengths = (0.01, 0.01, 0.01)

        #----default values for positions and separations of the grid-------------
        self.coords = [(0.005,) for _ in range(3)] # Coordenadas de los centros de los volúmenes
        self.dominios = [tuple([0.]+[self.coords[i][0]]+[self.lengths[i]]) for i in range(3)]
        self.deltas = [(0,) for _ in range(3)]
        # Hasta aquí tenemos un cubito

        self.__tags = {} # El etiquetado de todos los nodos sin fronteras
        self.__tags_fronteras = {} # El etiquetado de las fronteras
        
        self.dim = dim
        #---if a parameter is not a tuple, that parameter is transformed into a tuple
        if isinstance(volumes, int):  self.volumes = (volumes, 1, 1)
        if isinstance(lengths, (int, float)):  self.lengths = (lengths, lengths/10, lengths/10)
        
        # Si los parámetros son tuplas (pero no necesariamente sería una tupla de 3), arreglamos eso:
        if isinstance(volumes, tuple):
            faltan = 3 - len(volumes)
            self.volumes = volumes + tuple([1 for i in range(faltan)])
        if isinstance(lengths, tuple): 
            faltan = 3 - len(lengths)
            self.lengths = lengths + tuple([lengths[0]/10 for i in range(faltan)])
        #---------------------------------------------------------------
        
        # if volumes and lengths are given, initialize values acording to dimension
        if (volumes and lengths):       
            self.uniform_grid()
            self.init_tags()
            self.init_tags_fronteras()
    
    
    def uniform_grid(self):
        l = np.array(self.lengths) # Para el manejo con numpy
        v = np.array(self.volumes)
        d = l/v # Separación entre todos los nodos de cada dimensión
        start = d/2 # La frontera "inicial" del arreglo
        stop = l-d/2 # La frontera "final" del arreglo
        self.coords = [tuple(np.linspace(strt, stp, vol)) for strt, stp, vol in list(zip(start, stop, v))] # Meshgrid posible
        dominios = [np.insert(arr,(0,len(arr)),[0, l[idx]]) for idx, arr in enumerate(self.coords)] # Coordenadas + fronteras
        # Separación entre los nodos (Aquí hay que ver cómo es cuando tenemos un grid de 2x1x1 ya cuando se haga el FVM
        self.deltas = [self.set_deltas(dom)  if len(self.set_deltas(dom)) != 0 else (dom[-1],) for dom in dominios]
        self.dominios = [tuple(dom) for dom in dominios]
        #self.faces = [tuple((np.array(coords[:-1]) + np.array(coords[1:]))/2) for coords in self.coords]
        self.faces = [(self.dominios[idx][0],) + tuple((np.array(coords[:-1]) + np.array(coords[1:]))/2) \
                          + (self.dominios[idx][-1],) for idx, coords in enumerate(self.coords)]
        self.get_deltas_faces()
        self.get_grids()
    
    
    def set_deltas(self, dominio):
        """
        Método para obtener la distancia que hay entre los nodos
        """
        return tuple((dominio[1:]-dominio[:-1])[1:-1])


    # Creo que esto no se usará, pero estuvo chida la deducción, lo dejo de todos modos xd
    #def totalDomNodes(self):
    #    d_1 = 6*self.volumes[0] + 1
    #    d_2 = self.volumes[1]*d_1 - self.volumes[0]*(self.volumes[1] - 1)
    #    d_3 = self.volumes[2]*d_2 - self.volumes[0]*self.volumes[1]*(self.volumes[2] - 1)
    #    return d_3
    
    def init_tags(self):
        """
        Método que etiqueta las caras adyacentes de cada volumen dependiendo de la geometría. Pone un {} 
        (diccionario vacío) cuando es una frontera, una 'F' cuando es una cara interna y un 'Off' cuando 
        no se está contando esa cara por las dimensiones del problema. 
        """
        self.__tags = {}
        X, Y, Z = self.volumes
        for z in range(1,Z+1):
            for y in range(1,Y+1):
                for x in range(1,X+1):                   
                    t = b = n = s = "Off"
                    e = w = "F"
                    if x == 1: w = {}
                    elif x == X: e = {}
                    if self.dim > 1:
                        n = s = "F"
                        if y == 1: s = {}
                        elif y == Y: n = {}
                        if self.dim == 3:
                            t = b = "F"
                            if z == 1: b = {}
                            elif z == Z: t = {}
                    self.__tags[f"{x}-{y}-{z}"] = {"E": e, "W": w, "N": n, "S": s, "T": t, "B": b, 
                                             "coord": [self.dominios[0][x], self.dominios[1][y], self.dominios[2][z]]}
    
    def init_tags_fronteras(self):
        """
        Método para etiquetar las fronteras dependiendo de la dimensión, sólo se les da la propiedad de existir o no
        existir.
        """
        self.__tags_fronteras = {}
        X, Y, Z = [len(dom) for dom in self.dominios]
        for z in range(Z):
            for y in range(Y):
                for x in range(X):
                    t = b = n = s = "Off"
                    e = w = "ON"
                    if self.dim > 1: 
                        n = s = "ON"
                        if self.dim == 3: t = b = "ON"
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
                            self.__tags_fronteras[f"{x}-{y}-{z}"] = {"frontera": {var: value},
                                                 "coord": [self.dominios[0][x], self.dominios[1][y], self.dominios[2][z]],
                                                                  "cond": {}} 
                        elif z != 0 and z != (Z - 1):
                            if x != 0 and x != (X - 1):
                                if y == 0: var = "S"; value = s
                                elif y == (Y - 1) : var = "N"; value = n
                                self.__tags_fronteras[f"{x}-{y}-{z}"] = {"frontera": {var: value},
                                                 "coord": [self.dominios[0][x], self.dominios[1][y], self.dominios[2][z]],
                                                                      "cond": {}} 
                        else: continue
            
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
                
                
    def tag_wall_insulated(self, direction, coords=None):
        """
        Método para etiquetar fronteras con condición de aislamiento.
        """
        if coords:
            for idx, key in enumerate(coords):
                if key in list(self.__tags.keys()):
                    self.__tags[key][direction[idx]]["I"] = None
                elif key in list(self.__tags_fronteras.keys()):
                    self.__tags_fronteras[key]["cond"]["I"] = None
        else:
            if isinstance(direction, list):
                for idx, direct in enumerate(direction):
                    self.tag_wall(direct, "I", None)
            else:
                self.tag_wall(direction, "I", None)
                
    
    def set_dominio(self, dominio, faces=None):
        """
        Método para definir el dominio de estudio dadas unas coordenadas en forma de tupla.
        """
        # Si 'dominio' no es tupla, transforma 'dominio' a la tupla unidimensional (dominio,)
        # Tendría que ser una tupla de tuplas/listas/arreglos para que sea válido.
        if not isinstance(dominio, (tuple, int, float)): # Creo que si es una lista o un arreglo, no funciona enteros o float
            tupla = (tuple(dominio), self.dominios[1], self.dominios[2])
            dominio = tupla
        # Asigna los atributos de la mesh correspondientes    
        self.dominios = [tuple(dominio[i]) for i in range(3)]
        self.coords = [tuple(dominio[i][1:-1]) for i in range(3)]
        self.lengths = tuple([dominio[i][-1] for i in range(3)])
        self.volumes = tuple([len(dominio[i][1:-1]) for i in range(3)])
        #self.deltas = [self.set_deltas(np.array(dominio[i])) for i in range(3)]
        self.deltas = [self.set_deltas(np.array(dom))  if len(self.set_deltas(np.array(dom))) != 0 else (dom[-1],) for dom in self.dominios]
        
        if faces: 
            # Si me está pasando una lista (o sea, es de una dimensión)
            if isinstance(faces[0], (int, float)):
                self.faces = tuple(dominio[0]) + tuple(faces) + tuple(dominio[-1])
            else: # Suponemos aquí que nos está pasando una lista de listas (o tupla de tuplas)
                for idx, face_1dim in enumerate(faces):
                    faces[idx] = (dominio[idx][0], ) + tuple(face_1dim) + (dominio[idx][-1], )
                self.faces = faces
        else: 
            self.faces = [(self.dominios[idx][0],) + tuple((np.array(coords[:-1]) + np.array(coords[1:]))/2) \
                          + (self.dominios[idx][-1],) for idx, coords in enumerate(self.coords)]
        self.get_deltas_faces()
        self.init_tags()
        self.init_tags_fronteras()
        self.get_grids()
    
    def get_grid_deltas_dominios(self, axis=0, orientation="E"):
        deltas_dominios = []
        self.grid_deltas_dominios = np.array([])
        dominios = [np.array(dom) for dom in self.dominios]
        coords = [np.array(coord) for coord in self.coords]
        for direction in range(3):
            if direction != axis:
                if len(dominios[direction]) == 3:
                    deltas_dominios.append(dominios[direction][1])
                else:
                    deltas_dominios.append(coords[direction])
            else:
                deltas_dominios.append(dominios[direction][1:] - dominios[direction][:-1])
        if orientation == "E" or orientation == "N" or orientation == "T":
            deltas_dominios[axis] = deltas_dominios[axis][1:]
        else:
            deltas_dominios[axis] = deltas_dominios[axis][:-1]
        self.grid_deltas_dominios = np.meshgrid(deltas_dominios[0], deltas_dominios[1], deltas_dominios[2], 
                                                indexing='ij')
        return self.grid_deltas_dominios
    
    
    def get_deltas_faces(self):
        self.deltas_faces = []
        faces = [np.array(caras) for caras in self.faces]
        dominio = [np.array(doms) for doms in self.dominios]
        for direction in range(3):
            self.deltas_faces.append(faces[direction][1:] - faces[direction][:-1])
            
    def get_grids(self):
        self.grid_deltas = np.array([])
        self.grid_faces = np.array([])
        self.grid_coords = np.array([])
        deltas = self.deltas
        coords = self.coords
        faces = self.faces
        self.grid_deltas = np.meshgrid(deltas[0], deltas[1], deltas[2], indexing='ij')
        self.grid_coords = np.meshgrid(coords[0], coords[1], coords[2], indexing='ij')
        self.grid_faces = np.meshgrid(faces[0], faces[1], faces[2], indexing='ij')
                
    
    def info(self):
        """
        Método para imprimir información relevante del mallado
        """
        print('=====================================')
        print('     MESH INFORMATION   ')
        print('=====================================')
        print("\nMesh type: Cartesian")
        print(f"Number of dimensions of mesh: {self.dim}")
        variables = "X Y Z".split()
        for idx in range(self.dim):
            var = variables[idx]
            print(f"\n ----- {var}-axis -----")
            print(f"Number of {var.lower()} volumes: {self.volumes[idx]}")
            print(f"Lenght {var.lower()} of problem domain: {self.lengths[idx]}")
            print(f"List of {var.lower()} positions of volumes: \n{self.coords[idx]}")
            print(f"List of {var.lower()} positions of domain nodes: \n{self.dominios[idx]}")
            
            
    def draw(self):
        """
        Método para graficar la malla. Este método se invoca hasta que se hayan inizializado todas las condiciones
        de frontera.
        """
        # Graficamos las fronteras, sean o no activas
        dic_colors = {"D": "darkturquoise", "N": "red", "S": "magenta", "Off": "white", "I": "gray"}
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


    def get_area(self, direction = 0):
        """
        Método que regresa las áreas del volumen en la dirección indicada
        """
        perpendicular = [i for i in range(3) if i != direction]
        num_fronteras = self.volumes[direction]
        arreglo = [np.array([]) for _ in range(3)]
        arreglo[direction] = np.ones(num_fronteras)
        for i in perpendicular:
            arreglo[i] = self.deltas_faces[i]
        areas_grid = np.meshgrid(arreglo[0], arreglo[1], arreglo[2], indexing='ij')
        
        return areas_grid[perpendicular[0]]*areas_grid[perpendicular[1]]
    
    def stag_def(self, direct=0, stag="b"):
        """Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the X axis. """
        strategy_stag = lambda x: [self.dominios[direct][0]]+list(x)[:-1]+[self.coords[direct][-1]] if stag=="b" \
                                    else [self.coords[direct][0]]+list(x)[1:]+[self.dominios[direct][-1]]
        stg = strategy_stag(self.faces[direct])
        arreglo = [stg if idx==direct else self.dominios[idx] for idx in range(3)]
        strategy_stag_faces = lambda x: self.coords[x][:-1] if stag=="b" else self.coords[x][1:]
        stg_faces = strategy_stag_faces(direct)
        arreglo2 = [stg_faces if idx==direct else self.coords[idx] for idx in range(3)]
        return (arreglo[0], arreglo[1], arreglo[2]), (arreglo2[0], arreglo2[1], arreglo2[2])