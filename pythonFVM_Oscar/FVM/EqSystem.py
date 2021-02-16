#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:55:07 2019

@author: jose
"""


import numpy as np
import scipy.sparse as sp
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix

class EqSystem():
    """
    Class that uses a Coefficient object to construct the matrix  A, and vector
    b, for the system of equations Ax=b that represents the discrete version
    of the physical problem; this EqSystem class has methods and attributes
    for that putpose.
    
    
    Methods:

       
    Attribures:

    """
    
    def __init__(self, coeff):
        self._coeff = coeff
        malla = coeff.mesh()
        self._dim = malla.dim()
        self._dirichTagDic = malla.dirichValues()
        self._neumTagDic = malla.neumValues()
        self._nvx = malla.nvx() ; self._nvy = malla.nvy() ; self._nvz = malla.nvz()
        self._N = self._nvx*self._nvy*self._nvz
        self._A = None ; self._b = None
        
    def NaiveSetMatrix(self):
        if self._dim == 1: self.NaiveSetMatrix1D()
        elif self._dim == 2: self.NaiveSetMatrix2D()
        elif self._dim == 3: self.NaiveSetMatrix3D()
        
    def NaiveSetCSR(self, numb=False):
        if numb:
            if self._dim==3: self.NaiveSetCSR3D_numb()
            elif self._dim==2: self.NaiveSetCSR2D_numb()
            elif self._dim==1: self.NaiveSetCSR1D_numb()
        else:
            if self._dim==3: self.NaiveSetCSR3D()
            elif self._dim==2: self.NaiveSetCSR2D()
            elif self._dim==1: self.NaiveSetCSR1D()
            
    def NaiveSetCSC(self):
        self.NaiveSetCSR()
        self._A = csc_matrix(self.mat())

        
    def setMatrix_CSR(self):
        if self._dim==3: self.setMatrix3D_CSR()
        elif self._dim==2: self.setMatrix2D_CSR()
        elif self._dim==1: self.setMatrix1D_CSR()
 
    def setEfMatrix_CSR(self):
        if self._dim==3: self.setEfMatrix3D_CSR()
        elif self._dim==2: self.setEfMatrix2D_CSR()
        elif self._dim==1: self.setEfMatrix1D_CSR()       
        

    def NaiveSetMatrix1D(self):
        
        """Builds the matrix 'A' and the vector 'b' for a 1D problem using the 
        'coeff' attribute. When filling the matrix, for every volume the method 
        checks whether the volume has a border and if so, makes the respective
        correction to the corresponding A and b index."""
        
        coeff = self._coeff
        malla = coeff.mesh()
        tags = malla.tags()
        N = coeff.N()
        b = coeff.Su() #vector b del sistema de ec's Ax=b
        b = b.flatten()
        self._A = np.eye(N)
        A = self._A
        aP = coeff.aP()
        aE = coeff.aE() ; aW = coeff.aW()
        dx = malla.deltaX() 
        
        dxEast = dx[1:] ; dxWest = dx[:-1]
        Nx = malla.NdomX() ; Ny = malla.NdomY()
        nvx = malla.nvx() ; nvy = malla.nvy()
        
        for mesh_index in malla.unkown_mesh_indx():
            (i_mesh, j_mesh, k_mesh) = malla.meshIndx_ijk(mesh_index)
            Wtag = tags[(i_mesh-1) + j_mesh*Nx + k_mesh*Nx*Ny]
            Etag = tags[(i_mesh+1) + j_mesh*Nx + k_mesh*Nx*Ny]
            (i_vol, j_vol, k_vol) = (i_mesh - 1, j_mesh - 1, k_mesh - 1)
            vol_index = i_vol + j_vol*nvx + k_vol*nvx*nvy
            #------------------ FRONTERAS X ----------------------
            if not((Wtag == 'I') | (Wtag[0] == 'S')): #si se esta en la frontera West 
                (correc_b,correc_aP) = self.correcConst1(Wtag, aW[k_vol, j_vol, i_vol], dxWest[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol, j_vol, i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index, vol_index-1] = aW[k_vol, j_vol, i_vol]
                
            if not((Etag == 'I') | (Etag[0] == 'S')):  #Si se esta en la frontera Eastne
                (correc_b, correc_aP) = self.correcConst2(Etag, aE[k_vol, j_vol, i_vol], dxEast[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol, j_vol, i_vol] += correc_aP
            else:
                A[vol_index, vol_index+1] = aE[k_vol, j_vol, i_vol]
            A[vol_index, vol_index]=aP[k_vol, j_vol, i_vol]
        self._b = b

    def NaiveSetMatrix2D(self):
        """Builds the matrix 'A' and the vector 'b' for a 2D problem using the 
        'coeff' attribute. When filling the matrix, for every volume the method 
        checks whether the volume has a border and if so, makes the respective
        correction to the corresponding A and b index."""        
        coeff = self._coeff
        malla = coeff.mesh()
        tags=malla.tags()
        N = coeff.N()
        b = coeff.Su() #vector b del sistema de ec's Ax=b
        b = b.flatten()
        self._A = np.eye(N)
        A = self._A
        aP = coeff.aP()
        aE = coeff.aE() ; aW = coeff.aW()
        aN = coeff.aN() ; aS = coeff.aS()
        dx = malla.deltaX() ; dy = malla.deltaY()
        
        dxEast = dx[1:] ; dxWest = dx[:-1]
        dyNorth = dy[1:] ; dySouth = dy[:-1]
        Nx=malla.NdomX() ; Ny=malla.NdomY()
        nvx=malla.nvx() ; nvy = malla.nvy()
        
        for mesh_index in malla.unkown_mesh_indx():
            
            (i_mesh,j_mesh,k_mesh) = malla.meshIndx_ijk(mesh_index)
            
            Wtag=tags[(i_mesh-1) + j_mesh*Nx + k_mesh*Nx*Ny]
            Etag=tags[(i_mesh+1) + j_mesh*Nx + k_mesh*Nx*Ny]
            Stag=tags[i_mesh + (j_mesh-1)*Nx + k_mesh*Nx*Ny]
            Ntag=tags[i_mesh + (j_mesh+1)*Nx + k_mesh*Nx*Ny]
            (i_vol , j_vol , k_vol) = (i_mesh-1 , j_mesh-1 , k_mesh-1)
            vol_index = i_vol + j_vol*nvx + k_vol*nvx*nvy
            #------------------ FRONTERAS X ----------------------
            if not((Wtag=='I') | (Wtag[0]=='S')): #si se esta en la frontera West 
                (correc_b,correc_aP) = self.correcConst1(Wtag,aW[k_vol,j_vol,i_vol],dxWest[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index,vol_index-1] = aW[k_vol,j_vol,i_vol]
            if not((Etag=='I') | (Etag[0]=='S')):  #Si se esta en la frontera Eastne
                (correc_b,correc_aP) = self.correcConst2(Etag,aE[k_vol,j_vol,i_vol],dxEast[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else:
                A[vol_index,vol_index+1] = aE[k_vol,j_vol,i_vol]
            #------------------- Fronteras Y ------------------------------------------
            if not((Stag=='I') | (Stag[0]=='S')): #si se esta en la frontera South
                (correc_b,correc_aP) = self.correcConst1(Stag,aS[k_vol,j_vol,i_vol],dySouth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index,vol_index-nvx] = aS[k_vol,j_vol,i_vol]
            if not((Ntag=='I') | (Ntag[0]=='S')):  #Si se esta en la frontera North
                (correc_b,correc_aP) = self.correcConst2(Ntag,aN[k_vol,j_vol,i_vol],dyNorth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else:
                A[vol_index,vol_index+nvx] = aN[k_vol,j_vol,i_vol]
            
            A[vol_index,vol_index]=aP[k_vol,j_vol,i_vol]
        self._b = b
        
    def NaiveSetMatrix3D(self):
        """Builds the matrix 'A' and the vector 'b' for a 1D problem using the 
        'coeff' attribute. When filling the matrix, for every volume the method 
        checks whether the volume has a border and if so, makes the respective
        correction to the corresponding A and b index."""        
        coeff = self._coeff
        malla = coeff.mesh()
        tags=malla.tags()
        N = coeff.N()
        b = coeff.Su() #vector b del sistema de ec's Ax=b
        b = b.flatten()
        self._A = np.eye(N)
        A = self._A
        aP = coeff.aP()
        aE = coeff.aE() ; aW = coeff.aW() ; aT=coeff.aT()
        aN = coeff.aN() ; aS = coeff.aS() ; aB=coeff.aB()
        dx = malla.deltaX() ; dy = malla.deltaY() ; dz= malla.deltaZ()
        
        dxEast = dx[1:] ; dxWest = dx[:-1]
        dyNorth = dy[1:] ; dySouth = dy[:-1]
        dzTop = dz[1:] ; dzBottom = dz[:-1]
        Nx=malla.NdomX() ; Ny=malla.NdomY() 
        nvx=malla.nvx() ; nvy = malla.nvy() 
        
        for mesh_index in malla.unkown_mesh_indx():
            
            (i_mesh,j_mesh,k_mesh) = malla.meshIndx_ijk(mesh_index)
            
            Wtag=tags[(i_mesh-1) + j_mesh*Nx + k_mesh*Nx*Ny]
            Etag=tags[(i_mesh+1) + j_mesh*Nx + k_mesh*Nx*Ny]
            Stag=tags[i_mesh + (j_mesh-1)*Nx + k_mesh*Nx*Ny]
            Ntag=tags[i_mesh + (j_mesh+1)*Nx + k_mesh*Nx*Ny]
            Btag=tags[i_mesh + j_mesh*Nx + (k_mesh-1)*Nx*Ny]
            Ttag=tags[i_mesh + j_mesh*Nx + (k_mesh+1)*Nx*Ny]
            (i_vol , j_vol , k_vol) = (i_mesh-1 , j_mesh-1 , k_mesh-1)
            vol_index = i_vol + j_vol*nvx + k_vol*nvx*nvy
            #------------------ FRONTERAS X ----------------------
            if not((Wtag=='I') | (Wtag[0]=='S')): #si se esta en la frontera West 
                (correc_b,correc_aP) = self.correcConst1(Wtag,aW[k_vol,j_vol,i_vol],dxWest[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index,vol_index-1] = aW[k_vol,j_vol,i_vol]
            if not((Etag=='I') | (Etag[0]=='S')):  #Si se esta en la frontera Eastne
                (correc_b,correc_aP) = self.correcConst2(Etag,aE[k_vol,j_vol,i_vol],dxEast[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else:
                A[vol_index,vol_index+1] = aE[k_vol,j_vol,i_vol]
            #------------------- Fronteras Y ------------------------------------------
            if not((Stag=='I') | (Stag[0]=='S')): #si se esta en la frontera South
                (correc_b,correc_aP) = self.correcConst1(Stag,aS[k_vol,j_vol,i_vol],dySouth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index,vol_index-nvx] = aS[k_vol,j_vol,i_vol]
            if not((Ntag=='I') | (Ntag[0]=='S')):  #Si se esta en la frontera North
                (correc_b,correc_aP) = self.correcConst2(Ntag,aN[k_vol,j_vol,i_vol],dyNorth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else:
                A[vol_index,vol_index+nvx] = aN[k_vol,j_vol,i_vol]
            #------------------- Fronteras Z ------------------------------------------
            if not((Btag=='I') | (Ttag[0]=='S')): #si se esta en la frontera South
                (correc_b,correc_aP) = self.correcConst1(Btag,aB[k_vol,j_vol,i_vol],dzBottom[k_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                A[vol_index,vol_index-(nvx*nvy)] = aB[k_vol,j_vol,i_vol]
            if not((Ttag=='I') | (Ttag[0]=='S')):  #Si se esta en la frontera North
                (correc_b,correc_aP) = self.correcConst2(Ttag,aT[k_vol,j_vol,i_vol],dzTop[k_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                b[vol_index] += correc_b
                aP[k_vol,j_vol,i_vol] += correc_aP
            else:
                A[vol_index,vol_index+(nvx*nvy)] = aT[k_vol,j_vol,i_vol]  
                
            A[vol_index,vol_index]=aP[k_vol,j_vol,i_vol]
        self._b = b

    def NaiveSetCSR2D(self):
        
        coeff = self._coeff
        malla = coeff.mesh()
        tags=malla.tags()
        N = coeff.N()
        b = coeff.Su() #vector b del sistema de ec's Ax=b
        b = b.flatten()
        nvz,nvy,nvx = coeff.nvzyx()
        num_aE = (nvx-1)*nvy*nvz
        num_aN = nvx*(nvy-1)*nvz
        num_aT = nvx*nvy*(nvz-1)
        non_zeros = np.zeros( N +2*(num_aE+num_aN+num_aT) )
        JA = np.zeros( N +2*(num_aE+num_aN+num_aT) )
        IA = np.zeros( N+1 )
        aP = coeff.aP()
        aE = coeff.aE() ; aW = coeff.aW()
        aN = coeff.aN() ; aS = coeff.aS()
        dx = malla.deltaX() ; dy = malla.deltaY()
        
        dxEast = dx[1:] ; dxWest = dx[:-1]
        dyNorth = dy[1:] ; dySouth = dy[:-1]
        Nx=malla.NdomX() ; Ny=malla.NdomY()
        
        nzero_count =0 ; IA_count =1
        for k_vol in range(nvz):
            for j_vol in range(nvy):
                for i_vol in range(nvx):
            
                    (i_mesh , j_mesh , k_mesh) = (i_vol+1 , j_vol+1 , k_vol+1)
                    
                    Wtag=tags[(i_mesh-1) + j_mesh*Nx + k_mesh*Nx*Ny]
                    Etag=tags[(i_mesh+1) + j_mesh*Nx + k_mesh*Nx*Ny]
                    Stag=tags[i_mesh + (j_mesh-1)*Nx + k_mesh*Nx*Ny]
                    Ntag=tags[i_mesh + (j_mesh+1)*Nx + k_mesh*Nx*Ny]
                    vol_index = i_vol + j_vol*nvx + k_vol*nvx*nvy
                    
                    tag=tags[i_mesh + j_mesh*Nx + k_mesh*Nx*Ny]
                    if ((tag=='I') | (tag=='S')):
                        
                        scalar_aP = aP[k_vol,j_vol,i_vol]
                        #------------------ FRONTERAS X ----------------------
                        if not((Wtag=='I') | (Wtag[0]=='S')): #si se esta en la frontera West 
                            (correc_b,correc_aP) = self.correcConst1(Wtag,aW[k_vol,j_vol,i_vol],dxWest[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                            non_zeros[nzero_count] = aW[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index-1
                            nzero_count +=1
                        if not((Etag=='I') | (Etag[0]=='S')):  #Si se esta en la frontera Eastne
                            (correc_b,correc_aP) = self.correcConst2(Etag,aE[k_vol,j_vol,i_vol],dxEast[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else:
                            non_zeros[nzero_count] = aE[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index+1
                            nzero_count +=1
                        #------------------- Fronteras Y ------------------------------------------
                        if not((Stag=='I') | (Stag[0]=='S')): #si se esta en la frontera South
                            (correc_b,correc_aP) = self.correcConst1(Stag,aS[k_vol,j_vol,i_vol],dySouth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                            non_zeros[nzero_count] = aS[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index-nvx
                            nzero_count +=1
                        if not((Ntag=='I') | (Ntag[0]=='S')):  #Si se esta en la frontera North
                            (correc_b,correc_aP) = self.correcConst2(Ntag,aN[k_vol,j_vol,i_vol],dyNorth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else:
                            non_zeros[nzero_count] = aN[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index+nvx
                            nzero_count +=1
                            
                        non_zeros[nzero_count] = scalar_aP
                        
                    else:
                        non_zeros[nzero_count] = 1.      
                    
                    JA[nzero_count]=vol_index
                    nzero_count += 1
                    IA[IA_count] = nzero_count
                    IA_count += 1
                   
        self._A = sp.csr_matrix( (non_zeros,JA,IA), shape=(N,N) )
        self._b = b

    def NaiveSetCSR3D(self):
        
        coeff = self._coeff
        malla = coeff.mesh()
        tags=malla.tags()
        N = coeff.N()
        b = coeff.Su() #vector b del sistema de ec's Ax=b
        b = b.flatten()
        nvz,nvy,nvx = coeff.nvzyx()
        num_aE = (nvx-1)*nvy*nvz
        num_aN = nvx*(nvy-1)*nvz
        num_aT = nvx*nvy*(nvz-1)
        non_zeros = np.zeros( N +2*(num_aE+num_aN+num_aT) )
        JA = np.zeros( N +2*(num_aE+num_aN+num_aT) )
        IA = np.zeros( N+1 )
        aP = coeff.aP()
        aE = coeff.aE() ; aW = coeff.aW()
        aN = coeff.aN() ; aS = coeff.aS()
        aT = coeff.aT() ; aB = coeff.aB()
        dx = malla.deltaX() ; dy = malla.deltaY() ; dz=malla.deltaZ()
        
        dxEast = dx[1:] ; dxWest = dx[:-1]
        dyNorth = dy[1:] ; dySouth = dy[:-1]
        dzTop = dz[1:] ; dzBottom = dz[:-1]
        Nx=malla.NdomX() ; Ny=malla.NdomY()
        
        nzero_count =0 ; IA_count =1
        for k_vol in range(nvz):
            for j_vol in range(nvy):
                for i_vol in range(nvx):
            
                    (i_mesh , j_mesh , k_mesh) = (i_vol+1 , j_vol+1 , k_vol+1)
                    
                    Wtag=tags[(i_mesh-1) + j_mesh*Nx + k_mesh*Nx*Ny]
                    Etag=tags[(i_mesh+1) + j_mesh*Nx + k_mesh*Nx*Ny]
                    Stag=tags[i_mesh + (j_mesh-1)*Nx + k_mesh*Nx*Ny]
                    Ntag=tags[i_mesh + (j_mesh+1)*Nx + k_mesh*Nx*Ny]
                    Btag=tags[i_mesh + j_mesh*Nx + (k_mesh-1)*Nx*Ny]
                    Ttag=tags[i_mesh + j_mesh*Nx + (k_mesh+1)*Nx*Ny]
                    vol_index = i_vol + j_vol*nvx + k_vol*nvx*nvy
                    
                    tag=tags[i_mesh + j_mesh*Nx + k_mesh*Nx*Ny]
                    if ((tag=='I') | (tag=='S')):
                        
                        scalar_aP = aP[k_vol,j_vol,i_vol]
                        #------------------ FRONTERAS X ----------------------
                        if not((Wtag=='I') | (Wtag[0]=='S')): #si se esta en la frontera West 
                            (correc_b,correc_aP) = self.correcConst1(Wtag,aW[k_vol,j_vol,i_vol],dxWest[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                            non_zeros[nzero_count] = aW[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index-1
                            nzero_count +=1
                        if not((Etag=='I') | (Etag[0]=='S')):  #Si se esta en la frontera Eastne
                            (correc_b,correc_aP) = self.correcConst2(Etag,aE[k_vol,j_vol,i_vol],dxEast[i_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else:
                            non_zeros[nzero_count] = aE[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index+1
                            nzero_count +=1
                        #------------------- Fronteras Y ------------------------------------------
                        if not((Stag=='I') | (Stag[0]=='S')): #si se esta en la frontera South
                            (correc_b,correc_aP) = self.correcConst1(Stag,aS[k_vol,j_vol,i_vol],dySouth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                            non_zeros[nzero_count] = aS[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index-nvx
                            nzero_count +=1
                        if not((Ntag=='I') | (Ntag[0]=='S')):  #Si se esta en la frontera North
                            (correc_b,correc_aP) = self.correcConst2(Ntag,aN[k_vol,j_vol,i_vol],dyNorth[j_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else:
                            non_zeros[nzero_count] = aN[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index+nvx
                            nzero_count +=1
                        #------------------- Fronteras Z ------------------------------------------
                        if not((Btag=='I') | (Btag[0]=='S')): #si se esta en la frontera South
                            (correc_b,correc_aP) = self.correcConst1(Btag,aB[k_vol,j_vol,i_vol],dzBottom[k_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else: # Si no es de frontera asigna a la columa correspondiente pero no realiza corrección
                            non_zeros[nzero_count] = aB[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index-(nvx*nvy)
                            nzero_count +=1
                        if not((Ttag=='I') | (Ttag[0]=='S')):  #Si se esta en la frontera North
                            (correc_b,correc_aP) = self.correcConst2(Ttag,aT[k_vol,j_vol,i_vol],dzTop[k_vol]) #Realiza las corrección, la correción 1 aplica para w,s,b mientras que la corrección 2 aplica para e,n,t
                            b[vol_index] += correc_b
                            scalar_aP += correc_aP
                        else:
                            non_zeros[nzero_count] = aT[k_vol,j_vol,i_vol]
                            JA[nzero_count]=vol_index+(nvx*nvy)
                            nzero_count +=1
            
                        non_zeros[nzero_count] = scalar_aP
                        
                    else:
                        non_zeros[nzero_count] = 1.      
                    
                    JA[nzero_count]=vol_index
                    nzero_count += 1
                    IA[IA_count] = nzero_count
                    IA_count += 1
                   
        self._A = sp.csr_matrix( (non_zeros,JA,IA), shape=(N,N) )
        self._b = b


        
    def correcConst1(self,tag,a,delta):
        
        neumTagDic=self._neumTagDic
        dirichTagDic=self._dirichTagDic
        tagType=tag[0]
        if tagType=='D':
            tagValue=dirichTagDic[tag]
            correc_aP=0
            correc_b=-a*tagValue
        elif tagType=='N': 
            tagValue=neumTagDic[tag]
            correc_aP=a
            correc_b=a*tagValue*delta            
        return correc_b,correc_aP

    def correcConst2(self,tag,a,delta):
        
        neumTagDic=self._neumTagDic
        dirichTagDic=self._dirichTagDic
        tagType=tag[0]
        if tagType=='D':
            tagValue=dirichTagDic[tag]
            correc_aP=0
            correc_b=-a*tagValue
        elif tagType=='N': 
            tagValue=neumTagDic[tag]
            correc_aP=a
            correc_b=-a*tagValue*delta            
        return correc_b,correc_aP    
    
    def mat(self):
        return self._A
    
    def b(self):
        return self._b