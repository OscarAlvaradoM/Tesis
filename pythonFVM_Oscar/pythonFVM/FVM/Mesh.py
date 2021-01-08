#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
@author: jose
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from itertools import product

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

        self.__dirichTagDict = {} # No importa por ahora
        self.__neumTagDict = {}   # No importa por ahora
        self.__sourcTagDict = {}  # No importa por ahora
        self.__tags = {} # El etiquetado de todos los nodos
        self.__tags_f = {} # El etiquetado de las caras entre volúmenes

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
            dom_nodes = self.totalDomNodes()
            self.init_tags()
            #Nx, Ny, Nz = [len(dom) for dom in self.__dominios]
#
            #self.__hasWallW =[1 + j*Nx + k*Nx*Ny for k, j in product(range(1, Nz-1), range(1, Ny-1))]
            #self.__hasWallE =[Nx-2 + j*Nx + k*Nx*Ny for k, j in product(range(1, Nz-1),range(1, Ny-1))]
            
            #if self.__dim > 1:
            #    self.__hasWallS =[i+(1)*Nx+k*Nx*Ny for k,i in product(range(1,Nz-1),range(1,Nx-1))]
            #    self.__hasWallN =[i+(Ny-2)*Nx+k*Nx*Ny for k,i in product(range(1,Nz-1),range(1,Nx-1))]
            #else:
            #    self.__hasWallS =[]
            #    self.__hasWallN =[]
                
            #if self.__dim > 2:
            #    self.__hasWallB =[i+j*Nx+(1)*Nx*Ny for j,i in product(range(1,Ny-1),range(1,Nx-1))]
            #    self.__hasWallT =[i+j*Nx+(Nz-2)*Nx*Ny for j,i in product(range(1,Ny-1),range(1,Nx-1))]
            #else:
            #    self.__hasWallB =[]
            #    self.__hasWallT =[]
    
    
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
        self.__faces = [tuple(np.array(self.__deltas[idx]) + self.__lengths[idx]) for idx, coord in enumerate(self.__coords)]
    
    def setDeltas(self, dominio):
        return tuple((dominio[1:]-dominio[:-1])[1:-1])
    
    
    def totalDomNodes(self):
        d_1 = 6*self.__volumes[0] + 1
        d_2 = self.__volumes[1]*d_1 - self.__volumes[0]*(self.__volumes[1] - 1)
        d_3 = self.__volumes[2]*d_2 - self.__volumes[0]*self.__volumes[1]*(self.__volumes[2] - 1)
        return d_3
    
    
    def init_tags(self):
        """
        Clase que etiqueta las caras adyacentes de cada volumen dependiendo de la geometría. Pone un 0 (cero) cuando es una 
        frontera, una 'F' cuando es una cara interna y un 'Off' cuando no se está contando esa cara por las dimensiones del 
        problema. 
        """
        X, Y, Z = self.__volumes
        for z in range(Z):
            for y in range(Y):
                for x in range(X):
                    t = b = n = s = "Off"
                    e = w = "F"
                    if x == 0: e = 0
                    elif x == X - 1: w = 0
                    if self.__dim > 1:
                        n = s = "F"
                        if y == 0: s = 0
                        elif y == Y - 1: n = 0
                        if self.__dim == 3:
                            t = b = "F"
                            if z == 0: b = 0
                            elif z == Z - 1: t = 0
                            
                    self.__tags[f"{x}{y}{z}"] = {"E": e, "W": w, "N": n, "S": s, "T": t, "B": b}
            
    
    def tagWestWall(self, tag, value):
        """
        This method defines the West Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the West Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float)
        """
        
        
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0] + 2
        Ny = self.__volumes[1] + 2
        Nz = self.__volumes[2] + 2
        tags = self.__tags
        i = 0
        for k in range(Nz):
            for j in range(Ny):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index] = tag
        
    
    def tagEastWall(self,tag,value):
        """This method defines the East Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the East Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float) """
                
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        Nz = self.__volumes[2]+2
        tags=self.__tags
        i=Nx-1
        for k in range(Nz):
            for j in range(Ny):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index]=tag

    def tagSouthWall(self,tag,value):
        """This method defines the South Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the South Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float) """
        
        
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        Nz = self.__volumes[2]+2
        tags=self.__tags
        j=0
        for k in range(Nz):
            for i in range(Nx):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index]=tag
                
    def tagNorthWall(self,tag,value):
        """This method defines the North Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the North Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float) """
        
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        Nz = self.__volumes[2]+2
        tags=self.__tags
        j=Ny-1
        for k in range(Nz):
            for i in range(Nx):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index]=tag
                
    def tagBottomWall(self,tag,value):
        """This method defines the Bottom Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the Bottom Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float) """
                
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        tags=self.__tags
        k=0
        for j in range(Ny):
            for i in range(Nx):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index]=tag

    def tagTopWall(self,tag,value):
        """This method defines the Top Wall border condition.
        This is done by getting the instance attibute 'tags' and setting it's elements 
        equal to the 'tag' given. In this process, just the elements corresponding
        to the Top Wall are modified.
        The tag must begin with a D,N or S as this represent a dirichlet, neumman or 
        source node; accordingly. For instance, a border that doesn't allow any flux
        of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string begiinning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding border condition (float) """
                
        self.allocValue(tag,value)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        Nz = self.__volumes[2]+2
        tags=self.__tags
        k=Nz-1
        for j in range(Ny):
            for i in range(Nx):
                global_index = i + j*Nx + k*Nx*Ny
                tags[global_index]=tag  
    
    
    def setDominio(self, dominio):
        
        self.__autoMesh = False # La posición de los nodos no se calcula en automático
        
        # Si 'dominio' no es tupla transforma 'dominio' a la tupla unidimensional (dominio,)
        # Tendría que ser una tupla de tuplas/listas/arreglos para que sea válido.
        if not isinstance(dominio, tuple):
            tupla = (dominio,)
            dominio = tupla
        # Asigna los atributos de la mesh correspondientes    
        self.__dominios = [tuple(dominio[i]) for i in range(self.__dim)]
        self.__coords = [tuple(dominio[i][1:-1]) for i in range(self.__dim)]
        self.__lengths = [tuple(dominio[i][-1]) for i in range(self.__dim)]
        self.__volumes = [len(dominio[i][1:-1]) for i in range(self.__dim)]
        self.__deltas = [self.setDeltas(dominio[i]) for i in range(self.__dim)]
            
        dom_nodes = self.totalDomNodes()
        self.__tags = ["I" for i in range(dom_nodes)]
        self.tagsWallsOff()
        
        #Nx, Ny, Nz = [len(dom) for dom in self.__dominios]
        #self.__hasWallW =[1+j*Nx+k*Nx*Ny for k,j in product(range(1,Nz-1),range(1,Ny-1))]
        #self.__hasWallE =[Nx-2+j*Nx+k*Nx*Ny for k,j in product(range(1,Nz-1),range(1,Ny-1))]
        #if self.__dim > 1:
        #    self.__hasWallS =[i+(1)*Nx+k*Nx*Ny for k,i in product(range(1,Nz-1),range(1,Nx-1))]
        #    self.__hasWallN =[i+(Ny-2)*Nx+k*Nx*Ny for k,i in product(range(1,Nz-1),range(1,Nx-1))]
        #else:
        #    self.__hasWallS =[]
        #    self.__hasWallN =[]
        
        #if self.__dim > 2:    
        #    self.__hasWallB =[i+j*Nx+(1)*Nx*Ny for j,i in product(range(1,Ny-1),range(1,Nx-1))]
        #    self.__hasWallT =[i+j*Nx+(Nz-2)*Nx*Ny for j,i in product(range(1,Ny-1),range(1,Nx-1))]
        #else:
        #    self.__hasWallB =[]
        #    self.__hasWallT =[]
    
    
    def info(self):
        #-----Prints relevant information about the mesh------
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
            

    def createCoord(self):
        """
        Create and stores the nodes coordinates
        """
        self.__emptyCoord = False # Sigo sin saber para qué sirve

        domX, domY, domZ = [dom for dom in self.__dominios]
        self.__coords_dom = [coord_dom.flatten() for coord_dom in np.meshgrid(domY, domZ, domX)]
        
        
    def delCoor(self):
        """
        Deletes the coordX, coordY and coordZ attributes
        """
        del(self.__coords_dom)
        
        
    def getCoord(self):
        """
        Get the nodes coordinates
        """
        if self.__emptyCoord: # Check if the coordinates are stored
            self.createCoords()
            coord = self.__coords_dom
        else:
            coord = self.__coords_dom
        return coord
        
        
    def coordList(self):
        """
        Get the coordinates in a printable fashion
        """
        coordX, coordY, coordZ = self.getCoord()
        coordList = np.transpose([coordX, coordY, coordZ]) # No estoy seguro si esto sirva, está raro
        return coordList
    
    def coordTags(self):
        """
        Get the nodes coordinates and it's respective tag
        """
        coordX, coordY, coordZ = self.getCoord()
        coordTags = (coordX, coordY, coordZ, self.__tags) # La notación de los tags está rara, hacer diccionario
        
        return coordTags
    
    def printCoordTags(self, includeOffs = False):
        """
        Method that gets the nodes coordinates and their respective tags and prints 
        them. For clarity, it won't print more than 400 nodes and the nodes with 'off'
        tags are not included by default. The off nodes can be included by setting
        includeOffs to True.
        
        incudeOffs: parameter used to show the nodes with 'off' tags (True/False)
        """
        
        (X, Y, Z, tags) = self.coordTags()
        Nx, Ny, Nz = [len(dom) for dom in self.__dominios]
        
        if len(X) < 400:
            print("The coordinates and tags are:")
            if includeOffs or self.__dim == 3:
                for i in range(len(X)): print(f"({X[i]:.3f},{Y[i]:.3f},{Z[i]:.3f}):{tags[i]}")
            
            else:
                if self.__dim ==1:
                    node_1st = Nx*Ny + Nx # First node that is used for that dimension
                    (X, tags) = (X[node_1st:node_1st + Nx], tags[node_1st:node_1st + Nx])
                    for i in range(len(X)): print(f"({X[i]:.3f}):{tags[i]}")
                elif self.__dim == 2:
                    node_1st = Nx*Ny #first node that is used for that dimension
                    (X, Y, tags) = (X[node_1st:node_1st + Nx*Ny], Y[node_1st:node_1st + Nx*Ny], 
                                    tags[node_1st:node_1st + Nx*Ny])
                    for i in range(len(X)): print(f"({X[i]:.3f},{Y[i]:.3f}):{tags[i]}")
        else: print("Too many nodes (>400) to print")
                
    
    def allocValue(self, tag, value):
        """
        Takes the 'tag'  and sets it as a key for the 'value' in the corresponding
        dictionary. If the tag begins with a D,N or S it's placed in the dirichlet, 
        neumman or source dictionary; accordingly. For instance, a border that doesn't allow
        any flux of the scalar propertie can be setted using a 'tag' like N1 or N5 or N30 or
        etc., with a value of zero, that is, a Neumann condition of 0.
        
        tag: string beginning with D, N or S to represent a dirichlet, neumann or source node (string).
        value: value of the corresponding boundary condition (float)
        """
        
        dirich = self.__dirichTagDict 
        neum = self.__neumTagDict
        sourc = self.__sourcTagDict
        if tag[0] == 'D': dirich[tag] = value
        elif tag[0] == 'N': neum[tag] = value
        elif tag[0] == 'S': sourc[tag] = value
                    
                    
    def setUsrTag(self,usrTag,usrValue,x_lim,y_lim=(1.,1.),z_lim=(1.,1.)):
        """ Sets a tag (usrTag) and it's value (usrValue) for all the nodes inside 
        a given rectangular prism. The rectangular prism is specified via the 
        minimum x, maximum x, minimum y, maximum y, minimum z, maximum z that the 
        prism intersects.
        
        maximum x - minimum x = length of the prism
        maximum y - minimum y = width of the prism
        maximum z - minimum z = height of the prism
        
        If the values of Z aren't given then the figure is a rectangle instead of a prism.
        If only the X values are given then the figure is a segment instead of a prism
        
        usrTag: Tag to be setted, must begin with D,N or S to represent dirichlet, neumann, or source node (string)
        usrValue: The corresponding value of the tag to be setted (float)
        x_lim: A tuple that contains minimum x and maximum x, in that order (tuple of floats)
        y_lim: A tuple that contains minimum y and maximum y, in that order (tuple of floats)
        z_lim: A tuple that contains minimum z and maximum z, in that order (tuple of floats)
        """
        self.allocValue(usrTag,usrValue)
        
        Nx = self.__volumes[0]+2
        Ny = self.__volumes[1]+2
        Nz = self.__volumes[2]+2
        tags=self.__tags
        
        x1, x2 = x_lim
        y1, y2 = y_lim
        z1, z2 = z_lim
        x1_index = np.argmax(x1<= self.__dominioX) # primer índice que cumple x>x1
        x2_index = Nx-np.argmax( np.flip(self.__dominioX<=x2) ) #último índice +1  que cumple x<x2
        y1_index = np.argmax(y1<= self.__dominioY)
        y2_index = Ny-np.argmax( np.flip(self.__dominioY<=y2) )
        z1_index = np.argmax(z1<= self.__dominioZ)
        z2_index = Nz-np.argmax( np.flip(self.__dominioZ<=z2) )

        for k in range(z1_index,z2_index):
            for j in range(y1_index,y2_index):
                for i in range(x1_index,x2_index):
                    global_index = i + j*Nx + k*Nx*Ny
                    tags[global_index]=usrTag
                    
    def tagByCond(self,usrTag,usrValue,cond_func):
        
        """ Sets a tag 'usrTag' and it's value 'usrValue' for all the nodes in 
        the domain that satisfies the condition definded by the function 'cond_func'.
        The 'cond_func' function must accept three parameters x,y,z and return a 
        True/False value. 
        
        f(x,y,z)  -----> True / False
        
        usrTag: Tag to be setted, must begin with D,N or S to represent dirichlet, neumann, or source node (string)
        usrValue: The corresponding value of the tag to be setted (float)
        cond_func: function to avaluate wheter the node positions satisfies the condition or not (function)
        """
        
        self.allocValue(usrTag,usrValue)
        
        Nx = self.__volumes[0]+2 ; Ny = self.__volumes[1]+2 ; Nz = self.__volumes[2]+2
        X = self.__dominioX ; Y = self.__dominioY ; Z = self.__dominioZ
        tags=self.__tags
        
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    x=X[i] ; y = Y[j] ; z= Z[k]
                    if cond_func(x,y,z):
                        global_index = i + j*Nx + k*Nx*Ny
                        tags[global_index]=usrTag
        
        
                    
    def unkown_mesh_indx(self):
        """Returns the indexes of all the nodes that doesn't have a boundary condition,
        that is, all the indexes corresponding to an I or S tag."""
        tags=self.__tags
        tags_serie=pd.Series(tags)
        bool_cond=tags_serie.str.contains('I|S', regex=True) #serie de pandas booleana que tiene True para los indices con I o S
        unkown_mesh_indx=tags_serie[bool_cond].index
        return unkown_mesh_indx
    
    
    def meshIndx_ijk(self,meshIndx):
        """Takes a given node index 'meshIndx', and returns the three-dimensionals
        triple i,j,k indexes. For instance the second node has index equal to one so 
        meshIndx_ijk(1) should return (1,0,0)
        
        meshIndx: index of a node in the defined mesh (int)
        """
        
        Nx=self.NdomX() ; Ny =self.NdomY()
        k_mesh=meshIndx//(Nx*Ny)
        layer_indx=meshIndx%(Nx*Ny)
        j_mesh=(layer_indx)//Nx
        i_mesh=(layer_indx)%Nx
        return (i_mesh,j_mesh,k_mesh)

    
    def fixWallVars(self,eq_sol,use_memry=True):
        """ """
        ### A change is needed for working with the 'S' tags
        #### to see how it could work see the extendSol() method
        dic = self.__neumTagDict.copy()
        dic.update(self.__dirichTagDict)
        tags = self.__tags
        intWallNodes = self.intWallNodes(use_memry)
        (intWallNodesI,intWallNodesJ,intWallNodesK) = self.meshIndx_ijk(intWallNodes)
        relevant_tags = [ tags[i] for i in intWallNodes ]
        relevant_values = [ dic[i] for i in relevant_tags]
        
        for tag,value,i,j,k in zip(relevant_tags,relevant_values,intWallNodesI,intWallNodesJ,intWallNodesK ):
            i_vol=i-1 ; j_vol=j-1 ; k_vol=k-1  # i ,j and k are mesh indexes so a -1 should be applied
            Nx=self.__volumes[0]+2 ; Ny=self.__volumes[1]+2 
            if tag[0]!='D':
                n=0
                flux=dic[tag]
                neighbour_tag=tag
                if tag[1]=='e':
                    while neighbour_tag!='I': 
                        n-=1
                        neighbour_tag=tags[ (i+n) + j*Nx + k*Nx*Ny  ]
                    delta=self.__X[i_vol]-self.__X[i_vol+n]
                    value=flux*delta + eq_sol[ k_vol , j_vol , i_vol+n ]
                elif tag[1]=='w':
                    while neighbour_tag!='I' :
                        n+=1
                        neighbour_tag=tags[ (i+n) + j*Nx + k*Nx*Ny  ]
                    delta=self.__X[i_vol]-self.__X[i_vol+n]
                    value=flux*delta + eq_sol[ k_vol , j_vol , i_vol+n ]
                elif tag[1]=='n':
                    while neighbour_tag!='I':
                        n-=1
                        neighbour_tag=tags[ i + (j+n)*Nx + k*Nx*Ny  ]
                    delta=self.__Y[j_vol]-self.__Y[j_vol+n]
                    value=flux*delta + eq_sol[ k_vol , j_vol+n , i_vol ]
                elif tag[1]=='s':
                    while neighbour_tag!='I':
                        n+=1
                        neighbour_tag=tags[ i + (j+n)*Nx + k*Nx*Ny  ]
                    delta=self.__Y[j_vol]-self.__Y[j_vol+n]
                    value=flux*delta + eq_sol[ k_vol , j_vol+n , i_vol ]
                elif tag[1]=='t':
                    while neighbour_tag!='I':
                        n-=1
                        neighbour_tag=tags[ i + j*Nx + (k+n)*Nx*Ny  ]
                    delta=self.__Z[k_vol]-self.__Z[k_vol+n]
                    value=flux*delta + eq_sol[ k_vol+n , j_vol , i_vol ]
                elif tag[1]=='b':
                    while neighbour_tag!='I':
                        n+=1
                        neighbour_tag=tags[ i + j*Nx + (k+n)*Nx*Ny  ]
                    delta=self.__Z[k_vol]-self.__Z[k_vol+n]
                    value=flux*delta + eq_sol[ k_vol+n , j_vol , i_vol ]
                
            eq_sol[k_vol,j_vol,i_vol]=value
        
        return eq_sol

    def extendSol(self,sol_int):
        """ """
        Nx=self.__volumes[0]+2 ; Ny=self.__volumes[1]+2 ; Nz=self.__volumes[2]+2
        T_extend=np.zeros( (Nz,Ny,Nx) )
        T_extend[1:-1,1:-1,1:-1] = sol_int
        dic = self.__neumTagDict.copy()
        dic.update(self.__dirichTagDict)
        tags = self.__tags
        extWallNodes = self.extWallNodes()
        (intWallNodesI,intWallNodesJ,intWallNodesK) = self.meshIndx_ijk(extWallNodes)
        relevant_tags = [ tags[i] for i in extWallNodes ]
        relevant_values = [ dic[i] for i in relevant_tags]
        
        for tag,value,i,j,k in zip(relevant_tags,relevant_values,intWallNodesI,intWallNodesJ,intWallNodesK ):
            #
            if tag[0]=='N':
                n=0
                flux=dic[tag]
                neighbour_tag=tag
                if tag[1]=='e':
                    while ((neighbour_tag!='I') & ((i+n)>0)): ##############
                        n-=1
                        neighbour_tag=tags[ (i+n) + j*Nx + k*Nx*Ny  ]
                    delta=self.__dominioX[i]-self.__dominioX[i+n]
                    value=flux*delta + T_extend[ k , j , i+n ]
                elif tag[1]=='w':
                    while ((neighbour_tag!='I') & ((i+n)<Nx)) :
                        n+=1
                        neighbour_tag=tags[ (i+n) + j*Nx + k*Nx*Ny  ]
                    delta=self.__dominioX[i]-self.__dominioX[i+n]
                    value=flux*delta + T_extend[ k , j , i+n ]
                elif tag[1]=='n':
                    while ((neighbour_tag!='I') & ((j+n)>0)):
                        n-=1
                        neighbour_tag=tags[ i + (j+n)*Nx + k*Nx*Ny  ]
                    delta=self.__dominioY[j]-self.__dominioY[j+n]
                    value=flux*delta + T_extend[ k , j+n , i ]
                elif tag[1]=='s':
                    while ((neighbour_tag!='I') & ((j+n)<Ny)):
                        n+=1
                        neighbour_tag=tags[ i + (j+n)*Nx + k*Nx*Ny  ]
                    delta=self.__dominioY[j]-self.__dominioY[j+n]
                    value=flux*delta + T_extend[ k , j+n , i ]
                elif tag[1]=='t':
                    while ((neighbour_tag!='I') & ((k+n)>0)):
                        n-=1
                        neighbour_tag=tags[ i + j*Nx + (k+n)*Nx*Ny  ]
                    delta=self.__dominioZ[k]-self.__dominioZ[k+n]
                    value=flux*delta + T_extend[ k+n , j , i ]
                elif tag[1]=='b':
                    while ((neighbour_tag!='I') & ((k+n)<Nz)):
                        n+=1
                        neighbour_tag=tags[ i + j*Nx + (k+n)*Nx*Ny  ]
                    delta=self.__dominioZ[k]-self.__dominioZ[k+n]
                    value=flux*delta + T_extend[ k+n , j , i ]
            elif tag[0]=='S':
                value = T_extend[k,j,i]
                
            T_extend[k,j,i]=value
        
        return T_extend

    def intWallNodes(self,use_memry=True):
        """ """
        nvx = self.__volumes[0] ; nvy = self.__volumes[1] ; nvz = self.__volumes[2]
        Nx=nvx+2 ; Ny=nvy+2 ; Nz=nvz+2 
        if use_memry:
            if self.__pressed==False:
                volsEqSys=[i+j*Nx+k*Nx*Ny for k,j,i in product(range(1,Nz-1),range(1,Ny-1), range(1,Nx-1)  ) ]
                tags=self.__tags
                tags_serie=pd.Series(tags)
                bool_cond=tags_serie.str.contains('D|N', regex=True) #serie de pandas booleana que tiene True para los indices con I o S
                wall_mesh_indx=tags_serie[bool_cond].index     
                conjunto1=set(volsEqSys)
                conjunto2=set(wall_mesh_indx)
                interseccion=conjunto1.intersection(conjunto2) #set
                intWallNodes=np.array( [element for element in interseccion] ) #array
                self.__intWallNodes = intWallNodes
            else:
                intWallNodes=self.__intWallNodes
            
        else:
            volsEqSys=[i+j*Nx+k*Nx*Ny for k,j,i in product(range(1,Nz-1),range(1,Ny-1), range(1,Nx-1)  ) ]
            tags=self.__tags
            tags_serie=pd.Series(tags)
            bool_cond=tags_serie.str.contains('D|N', regex=True) #serie de pandas booleana que tiene True para los indices con I o S
            wall_mesh_indx=tags_serie[bool_cond].index     
            conjunto1=set(volsEqSys)
            conjunto2=set(wall_mesh_indx)
            interseccion=conjunto1.intersection(conjunto2) #set
            intWallNodes=np.array( [element for element in interseccion] ) #array
   
        return intWallNodes

    def extWallNodes(self):
        """Get all indexes of the nodes that corresponds to the exterior walls,
        that is, the nodes of the W, E, N, S, N, T and B walls, according to
        the number of dimensions in the mesh"""
        nvx = self.__volumes[0] ; nvy = self.__volumes[1] ; nvz = self.__volumes[2]
        Nx=nvx+2 ; Ny=nvy+2 ; Nz=nvz+2 
        indxs2 = [] ; indxs3=[]
        if self.__dim ==1:
            indxs1= [Nx*(Ny+1), Nx*(Ny+2)-1]
        elif self.__dim==2:
            indxs1 = [i+j*Nx+ 1*Nx*Ny for j,i in product(range(1,Ny-1), [0,Nx-1] ) ]
            indxs2 = [i+j*Nx+ 1*Nx*Ny for j,i in product([0,Ny-1], range(0,Nx) ) ]
        elif self.__dim==3:
            indxs1 = [i+j*Nx+ k*Nx*Ny for k,j,i in product( [0,Nz-1],range(0,Ny), range(0,Nx) ) ]
            indxs2 = [i+j*Nx+ k*Nx*Ny for k,j,i in product(range(1,Nz-1),[0,Ny-1], range(0,Nx)  ) ]
            indxs3 = [i+j*Nx+k*Nx*Ny for k,j,i in product(range(1,Nz-1),range(1,Ny-1), [0,Nx-1]  ) ]
        indxs=np.array(indxs1+indxs2+indxs3)
        extWallNodes = np.sort(indxs)
        
        return extWallNodes
        
        
    def press(self):
        """This method states that the mesh won't require any changes and therefore
        some attributes will be saved in it's final value."""
        self.__pressed=True
    
    def uStagDef(self):
        """Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the X axis. """
        stgX=np.zeros(self.__volumes[0]+1)
        stgX[-1]=self.__lengths[0]
        stgX[1:-1]=0.5*(self.X()[:-1]+self.X()[1:])
        return (stgX,self.dominioY(),self.dominioZ())
    
    def vStagDef(self):
        """Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the Y axis. """
        
        stgY=np.zeros(self.__volumes[1]+1)
        stgY[-1]=self.__lengths[1]
        stgY[1:-1]=0.5*(self.Y()[:-1]+self.Y()[1:])
        return (self.dominioX(),stgY,self.dominioZ())
    
    def wStagDef(self):
        """Gives the parameters needed for defining a staggered mesh via the 'setDominio()'
        method. The staggering it's made in the Z axis. """
        
        stgZ=np.zeros(self.__volumes[2]+1)
        stgZ[-1]=self.__lengths[2]
        stgZ[1:-1]=0.5*(self.Z()[:-1]+self.Z()[1:])
        return (self.dominioX(),self.dominioY(),stgZ)
    
    def draw(self,marcador=None,brd_mrkr='s',new_fig=True, s_int=15, s_wall=50):
        """ Makes a plot of the mesh, that is, a plot which shows all the nodes positions.
        The dirichet, neumann, source and interior ones are represented with yellow, 
        blue, magenta and black; accordingly. In the 3D meshes this method is useful
        with meshes that have just a few nodes, for 3D meshes with many nodes the
        'drawByCut' method is preferred
        
        marcador: Marker for the non-border nodes (string)
        brd_mrkr: Marker for the border nodes (string)
        new_fig: Flag to decide if a new figure should be used (True/False)
        """
        
        (X,Y,Z,tags)=self.coordTags()
        Nx=self.NdomX() ; Ny=self.NdomY()
        if self.__dim ==1:
            node_1st = Nx*Ny+Nx #first node that is used for that dimension
            (X,tags)=(X[node_1st:node_1st+Nx],tags[node_1st:node_1st+Nx])
            data = {'X':X,'Tag':tags}
        elif self.__dim==2:
            node_1st = Nx*Ny #first node that is used for that dimension
            (X,Y,tags)=(X[node_1st:node_1st+Nx*Ny],Y[node_1st:node_1st+Nx*Ny],tags[node_1st:node_1st+Nx*Ny])
            data = {'X':X,'Y':Y,'Tag':tags}
        elif self.__dim ==3:
            data = {'X':X,'Y':Y,'Z':Z,'Tag':tags}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I', regex=False)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        inter_coordTags=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        bool_cond=df['Tag'].str.contains('D', regex=False)
        dirich_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('N', regex=False)
        neum_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('S', regex=False)
        source_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------
        if self.__dim ==1:
            if new_fig: plt.figure()
            plt.title("Grid")
            plt.xlim(0,self.__lengths[0])
            plt.scatter(dirich_coordTags[0],0*dirich_coordTags[0],color='y',marker=brd_mrkr,s=s_wall)
            plt.scatter(neum_coordTags[0],0*neum_coordTags[0],color='b',marker=brd_mrkr,s=s_wall)
            plt.scatter(source_coordTags[0],0*source_coordTags[0],marker=marcador,color='m',s=s_int)
            plt.scatter(inter_coordTags[0],0*inter_coordTags[0],marker=marcador,color='k',s=s_int)                
        elif self.__dim ==2:
            if new_fig: plt.figure()
            plt.title("Grid")
            plt.xlim(0,self.__lengths[0])
            plt.ylim(0,self.__lengths[1])
            plt.scatter(neum_coordTags[0],neum_coordTags[1],color='b',marker=brd_mrkr,s=s_wall)
            plt.scatter(dirich_coordTags[0],dirich_coordTags[1],color='y',marker=brd_mrkr,s=s_wall)
            plt.scatter(source_coordTags[0],source_coordTags[1],marker=marcador,color='m',s=s_int)
            plt.scatter(inter_coordTags[0],inter_coordTags[1],marker=marcador,color='k',s=s_int)

        elif self.__dim ==3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(dirich_coordTags[0], dirich_coordTags[1], dirich_coordTags[2], marker=brd_mrkr, color='y')
            ax.scatter(neum_coordTags[0],neum_coordTags[1] ,neum_coordTags[2] , marker=brd_mrkr,color='b')
            ax.scatter(source_coordTags[0],source_coordTags[1] , source_coordTags[2], marker=marcador,color='m',s=s_int)
            ax.scatter(inter_coordTags[0],inter_coordTags[1] , inter_coordTags[2], marker=marcador,color='k',s=s_int)
        plt.show()


    def drawByCut(self,x_cut=None,y_cut=None,z_cut=None):
        """ Makes three 2D plots of a 3D mesh, that is, three plots which shows the nodes
        positions corresponding to a transversal cutting plane. The three plots can be
        interpreted as some frontal, side and upper view.
        The x, y and z numbers that determine the location of the frontal, side and upper
        cutting plane are given by 'x_cut', 'y_cut' and 'z_cut' parameters, accordingly.
        The dirichet, neumann, source and interior ones are represented with yellow, 
        blue, magenta and black. 
        
        x_cut: x constant coordinate of the frontal cutting plane (float)
        y_cut: y constant coordinate of the side cutting plane (float)
        z_cut: z constant coordinate of the upper cutting plane (float)
        
        """

        Nx= self.NdomX() ; Ny = self.NdomY()  ; Nz = self.NdomZ()
        tags=self.tags() # Extraemmos los tags del dominio
        tags = np.reshape(tags,(Nz,Ny,Nx))
        lx = self.lx() ; ly = self.ly() ; lz = self.lz() #extraemos medidas del dominio
        domX = self.dominioX() ; domY = self.dominioY() ; domZ = self.dominioZ()
        coordY,coordZ,coordX = np.meshgrid(domY,domZ,domX)
        if x_cut==None: x_cut=lx/2.
        if y_cut==None: y_cut=ly/2.
        if z_cut==None: z_cut=lz/2.
        indx_xcut=np.argmax(x_cut<domX)
        indx_ycut=np.argmax(y_cut<domY)
        indx_zcut=np.argmax(z_cut<domZ)
        if (x_cut-domX[indx_xcut-1]) < (domX[indx_xcut] - x_cut): 
            indx_xcut-=1  #Si está mas cerca de X[índice anterior] tomamos el índice anterior
        if (y_cut-domY[indx_ycut-1]) < (domY[indx_ycut] - y_cut): 
            indx_ycut-=1
        if (z_cut-domZ[indx_zcut-1]) < (domZ[indx_zcut] - z_cut): 
            indx_zcut-=1
        #========================================================
        # GRAFICAMOS EL CORTE Y
        #========================================================
        X=coordX[:,indx_ycut,:] ; Y=coordY[:,indx_ycut,:] ; Z=coordZ[:,indx_ycut,:] 
        tagsCut = tags[:,indx_ycut,:]
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tagsCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I', regex=False)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        inter_coordTags=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        bool_cond=df['Tag'].str.contains('D', regex=False)
        dirich_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('N', regex=False)
        neum_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('S', regex=False)
        source_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------              
        plt.figure()
        plt.title("GridXZ   y={}".format(domY[indx_ycut]))
        plt.xlim(0,self.__lengths[0])
        plt.ylim(0,self.__lengths[2])
        plt.scatter(dirich_coordTags[0],dirich_coordTags[2],color='y',marker='s',s=50)
        plt.scatter(neum_coordTags[0],neum_coordTags[2],color='b',marker='s',s=50)
        plt.scatter(source_coordTags[0],source_coordTags[2],color='m',s=15)
        plt.scatter(inter_coordTags[0],inter_coordTags[2],color='k',s=15)
        plt.show()
        #========================================================
        # GRAFICAMOS EL CORTE X
        #========================================================
        X=coordX[:,:,indx_xcut] ; Y=coordY[:,:,indx_xcut] ;Z=coordZ[:,:,indx_xcut]
        tagsCut = tags[:,:,indx_xcut]
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tagsCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I', regex=False)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        inter_coordTags=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        bool_cond=df['Tag'].str.contains('D', regex=False)
        dirich_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('N', regex=False)
        neum_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('S', regex=False)
        source_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------              
        plt.figure()
        plt.title("GridYZ   x={}".format(domX[indx_xcut]))
        plt.xlim(0,self.__lengths[1])
        plt.ylim(0,self.__lengths[2])
        plt.scatter(dirich_coordTags[1],dirich_coordTags[2],color='y',marker='s',s=50)
        plt.scatter(neum_coordTags[1],neum_coordTags[2],color='b',marker='s',s=50)
        plt.scatter(source_coordTags[1],source_coordTags[2],color='m',s=15)
        plt.scatter(inter_coordTags[1],inter_coordTags[2],color='k',s=15)
        plt.show()
        #========================================================
        # GRAFICAMOS EL CORTE Z
        #========================================================
        X=coordX[indx_zcut,:,:] ; Y=coordY[indx_zcut,:,:] ;Z=coordZ[indx_zcut,:,:]
        tagsCut = tags[indx_zcut,:,:]
        data = {'X':X.flatten(),'Y':Y.flatten(),'Z':Z.flatten(),'Tag':tagsCut.flatten()}
        #-------Busquedas y clasificación-------------------------
        df = pd.DataFrame(data)
        bool_cond=df['Tag'].str.contains('I', regex=False)  #crea una serie de pandas booleana (a partir de la serie 'Tag') que es true en aquellos indices que contienen I
        inter_coordTags=np.transpose(df[bool_cond].values) #nos quedamos con un arreglo que tiene las coordenadas y tags de todos los puntos interiores
        bool_cond=df['Tag'].str.contains('D', regex=False)
        dirich_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('N', regex=False)
        neum_coordTags=np.transpose(df[bool_cond].values)
        bool_cond=df['Tag'].str.contains('S', regex=False)
        source_coordTags=np.transpose(df[bool_cond].values)
        #------------Graficación----------------------------------              
        plt.figure()
        plt.title("GridXY z={}".format(domZ[indx_zcut]))
        plt.xlim(0,self.__lengths[0])
        plt.ylim(0,self.__lengths[1])
        plt.scatter(dirich_coordTags[0],dirich_coordTags[1],color='y',marker='s',s=50)
        plt.scatter(neum_coordTags[0],neum_coordTags[1],color='b',marker='s',s=50)
        plt.scatter(source_coordTags[0],source_coordTags[1],color='m',s=15)
        plt.scatter(inter_coordTags[0],inter_coordTags[1],color='k',s=15)
        plt.show()


    def vols(self):
        """ Returns a 3-dimensional numpy array containing all the value of the volumes,
         that is lenght times width times height, for each one."""
        dim = self.__dim
        dx=self.deltaX() ; dy=self.deltaY() ; dz=self.deltaZ()
        dx_vol = (dx[1:]+dx[:-1])/2. 
        dy_vol = (dy[1:]+dy[:-1])/2.
        dz_vol = (dz[1:]+dz[:-1])/2.
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        vols_dim_grid = np.meshgrid(dy_vol,dz_vol,dx_vol)
        return vols_dim_grid[2]*vols_dim_grid[0]*vols_dim_grid[1]
    
    def areasX(self):
        """" Returns a 3-dimensional numpy array containing the areas of the faces in
        the X direction. The numpy array contains an area for each volume. """
        
        dim = self.__dim
        dy=self.deltaY() ; dz=self.deltaZ()
        nvx=self.nvx()
        dy_vol = (dy[1:]+dy[:-1])/2.  #longitud y (ancho) del volumen de control
        dz_vol = (dz[1:]+dz[:-1])/2.  #longitud z (alto) del volumen de control
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        areas_grid =np.meshgrid(dy_vol,dz_vol,np.ones(nvx) )
        return areas_grid[0]*areas_grid[1]

    def areasY(self):
        """" Returns a 3-dimensional numpy array containing the areas of the faces in
        the Y direction. The numpy array contains an area for each volume. """        

        dim = self.__dim
        dx=self.deltaX() ; dz=self.deltaZ()
        nvy=self.nvy()
        dx_vol = (dx[1:]+dx[:-1])/2.  #longitud y (ancho) del volumen de control
        dz_vol = (dz[1:]+dz[:-1])/2.  #longitud z (alto) del volumen de control
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >2: dz_vol[-1] += 0.5*dz[-1] ; dz_vol[0] += 0.5*dz[0]
        areas_grid =np.meshgrid(np.ones(nvy),dz_vol, dx_vol)
        return areas_grid[2]*areas_grid[1]

    def areasZ(self):
        """" Returns a 3-dimensional numpy array containing the areas of the faces in
        the Z direction. The numpy array contains an area for each volume. """  
        
        dim = self.__dim
        dx=self.deltaX() ; dy=self.deltaY()
        nvz=self.nvz()
        dx_vol = (dx[1:]+dx[:-1])/2.  #longitud y (ancho) del volumen de control
        dy_vol = (dy[1:]+dy[:-1])/2.  #longitud y (ancho) del volumen de control
        dx_vol[-1] += 0.5*dx[-1] ; dx_vol[0] += 0.5*dx[0]
        if dim >1: dy_vol[-1] += 0.5*dy[-1] ; dy_vol[0] += 0.5*dy[0]
        areas_grid =np.meshgrid(dy_vol,np.ones(nvz),dx_vol )
        return areas_grid[2]*areas_grid[0]

    def dxe(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's east neighbor. The numpy array contains 
        a separation (distance) for each volume. """  
        
        dx=self.deltaX()
        nvy=self.__volumes[1] ; nvz=self.__volumes[2]
        dxe = dx[1:]
        dxe_grid = np.meshgrid(np.ones(nvy),np.ones(nvz),dxe)
        return dxe_grid[2]

    def dxw(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's west neighbor. The numpy array contains 
        a separation (distance) for each volume. """    
        
        dx=self.deltaX()
        nvy=self.__volumes[1] ; nvz=self.__volumes[2]
        dxw = dx[:-1]
        dxw_grid = np.meshgrid(np.ones(nvy),np.ones(nvz),dxw)
        return dxw_grid[2]  

    def dyn(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's north neighbor. The numpy array contains 
        a separation (distance) for each volume. """     
        
        dy=self.deltaY()
        nvx=self.__volumes[0] ; nvz=self.__volumes[2]
        dyn = dy[1:]
        dyn_grid = np.meshgrid(dyn,np.ones(nvz),np.ones(nvx))
        return dyn_grid[0] 

    def dys(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's south neighbor. The numpy array contains 
        a separation (distance) for each volume. """     
        
        dy=self.deltaY()
        nvx=self.__volumes[0] ; nvz=self.__volumes[2]
        dys = dy[:-1]
        dys_grid = np.meshgrid(dys,np.ones(nvz),np.ones(nvx))
        return dys_grid[0]

    def dzt(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's top (upper) neighbor. The numpy array contains 
        a separation (distance) for each volume. """     
        
        dz=self.deltaZ()
        nvx=self.__volumes[0] ; nvy=self.__volumes[1]
        dzt = dz[1:]
        dzt_grid = np.meshgrid(np.ones(nvy),dzt,np.ones(nvx))
        return dzt_grid[1]

    def dzb(self):
        """" Returns a 3-dimensional numpy array containing the separations between 
        the corresponding node and it's bottom (lower) neighbor. The numpy array contains 
        a separation (distance) for each volume. """     
        
        dz=self.deltaZ()
        nvx=self.__volumes[0] ; nvy=self.__volumes[1]
        dzb = dz[:-1]
        dzb_grid = np.meshgrid(np.ones(nvy),dzb,np.ones(nvx))
        return dzb_grid[1]

    
    def xe(self,x=None,dxe=None):
        #Gives a 3-dimensional array containing the x coordinate of the east border in the volume.
        #The numpy array has a x coordinate (position) for every volume."""
        return x+0.5*dxe

    def xw(self,x=None,dxw=None):
        #Gives a 3-dimensional array containing the x coordinate of the west border in the volume.
        #The numpy array has a x coordinate (position) for every volume."""
        return x-0.5*dxw
    
    def yn(self,y=None,dyn=None):
        #Gives a 3-dimensional array containing the y coordinate of the north border in the volume.
        #The numpy array has a y coordinate (position) for every volume.
        return y+0.5*dyn

    def ys(self,y=None,dys=None):
        #Gives a 3-dimensional array containing the y coordinate of the south border in the volume.
        #The numpy array has a y coordinate (position) for every volume.
        return y-0.5*dys    

    def zt(self,z=None,dzt=None):
        #Gives a 3-dimensional array containing the z coordinate of the top border in the volume.
        #The numpy array has a z coordinate (position) for every volume.
        return z+0.5*dzt

    def zb(self,z=None,dzb=None):
        #Gives a 3-dimensional array containing the z coordinate of the bottom border in the volume.
        #The numpy array has a z coordinate (position) for every volume.
        return z-0.5*dzb        

    def nvx(self):
        return self.__volumes[0]

    def nvy(self):
        return self.__volumes[1]

    def nvz(self):
        return self.__volumes[2]
    
    def nvzyx(self):
        return (self.__volumes[2] , self.__volumes[1] , self.__volumes[0])

    def lx(self):
        return self.__lengths[0]

    def ly(self):
        return self.__lengths[1]

    def lz(self):
        return self.__lengths[2]
    
    def X(self):
        return self.__X

    def Y(self):
        return self.__Y

    def Z(self):
        return self.__Z

    def dominioX(self):
        return self.__dominioX

    def dominioY(self):
        return self.__dominioY

    def dominioZ(self):
        return self.__dominioZ

    def deltaX(self):
        return np.array(self.__deltaX)

    def deltaY(self):
        return np.array(self.__deltaY)

    def deltaZ(self):
        return np.array(self.__deltaZ)
    
    def dim(self):
        return self.__dim

    def coordGrid(self):
        coordGrid=np.meshgrid(self.Y(),self.Z(),self.X())
        return (coordGrid[2],coordGrid[0],coordGrid[1])
    
    def coordX(self):
        return self.__coordX

    def coordY(self):
        return self.__coordY

    def coordZ(self):
        return self.__coordZ
    
    def tags(self):
        return self.__tags
    
    def tagShaped(self):
        """Get the tags and their respective nodes coordinates in a printable fashion"""
        nvx=self.__volumes[0] ; nvy=self.__volumes[1] ; nvz=self.__volumes[2] 
        tags_shaped= np.reshape(self.__tags,(nvz,nvy,nvx)) #nvx+2
        return tags_shaped
    
    def dirichValues(self):
        return self.__dirichTagDict

    def neumValues(self):
        return self.__neumTagDict

    def sourcValues(self):
        return self.__sourcTagDict    
    
    

    
