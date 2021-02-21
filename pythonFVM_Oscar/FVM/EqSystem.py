from scipy.sparse import diags
import numpy as np

class EqSystem():
    """
    Class that uses a Coefficient object to construct the matrix  A, and vector
    b, for the system of equations Ax=b that represents the discrete version
    of the physical problem; this EqSystem class has methods and attributes
    for that putpose.
    
    
    Methods:

       
    Attribures:

    """
    
    def __init__(self, coef):
        self.coef = coef
        self.aP = coef.aP
        self.aE, self.aW = coef.aE, coef.aW
        self.aN, self.aS = coef.aN, coef.aS
        self.aT, self.aB = coef.aT, coef.aB
        self.Su, self.Sp = coef.Su, coef.Sp
        self.N = self.aP.shape[0]*self.aP.shape[1]*self.aP.shape[2]
        self.A, self.b = None, None
        
    def get_diag(self, array, k=0):
        """
        Método para construir una matriz diagonal dado un arreglo 2D y la diagonal en que queremos 
        poner ese arreglo
        """
        if self.N > 50:
            return diags(array, k)
        else:
            return np.diag(array,k)
        
        
    def get_A_matrix(self):
        """
        Método para construir la matriz A, de la ecuación Ax=b, a partir de los coeficientes obtenidos.
        """
        aP = np.ravel(self.aP)
        self.A = self.get_diag(aP)
        if self.coef.dim == 1:
            aE = np.ravel(self.aE)[:-1]
            aW = np.ravel(self.aW)[1:]
            self.A += self.get_diag(aE, k=1) + self.get_diag(aW, k=-1)
        elif self.coef.dim == 2:
            aN = np.ravel(self.aN)[:-1]
            aS = np.ravel(self.aS)[1:]
            aE = np.ravel(self.aE)[:-4]
            aW = np.ravel(self.aW)[4:]
            self.A += self.get_diag(aN, k=1) + self.get_diag(aS, k=-1) + self.get_diag(aE, k=4) + self.get_diag(aW, k=-4)
        # Esta parte está por verse, hay que encontrar un ejemplo chido en 3D
        if self.coef.dim == 3:
            pass
        return self.A
    
    
    def get_b_vector(self):
        """
        Método para obtener el vector b, de la ecuación Ax = b, a partir del arreglo Su construido anteriormente
        """
        self.b = np.ravel(self.Su)
        return self.b
    
    
    def get_solution(self):
        """
        Método para obtener los valores de la solución al sistema
        """
        A = self.get_A_matrix()
        b = self.get_b_vector()
        return np.linalg.solve(A, b)
