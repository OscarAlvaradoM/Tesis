#!/usr/bin/env python
# coding: utf-8

# # Notebook de benchmarking de la Ecuación de Poisson

# In[1]:


import numpy as np

import time
from csv import writer

import FVM as fvm
from FVM import Mesh


# ## Ecuación de Poisson 3D

# ### Construcción de la malla

# In[2]:


def benchmark_mesh(vol):
    mesh = fvm.Mesh(3, volumes = (vol, vol, vol), lengths = (vol/10, vol/10, vol/10))
    mesh.tag_wall_dirichlet("S", 100)
    mesh.tag_wall_dirichlet(["W", "E", "T", "N", "B"], [0,0,0,0,0])
    return mesh


# ### Haciendo cuentitas del FVM

# In[3]:


def benchmark_set_boundary_conditions(mesh):
    coef = fvm.Coefficients(mesh)
    coef.set_diffusion(1_000)
    return coef


# ### Obteniendo soluciones

# In[4]:


def benchmark_solutions(coef):
    sistema = fvm.EqSystem(coef)
    sistema.get_solution()


# ---

# ## Midiendo tiempo más formalmente

# In[5]:


def iterate_volumes(volumes, times):
    for volume in volumes:
        mesh = benchmark_mesh(volume)
        coefficients = benchmark_set_boundary_conditions(mesh)
        functions = [benchmark_mesh, benchmark_set_boundary_conditions, benchmark_solutions]
        args = [volume, mesh, coefficients]
        
        for f, arg in zip(functions, args):
            print(f"Comencé el de {volume} volúmenes con la función {f.__name__}")
            list_of_time_statistics = list_of_statistics(f, arg, volume, times)
            write_to_file(f, list_of_time_statistics)
        print(f"Terminé el de {volume} volúmenes")


# In[6]:


def list_of_statistics(f, arg, volume, times):
    time_list, μ, σ = get_statistics(f, arg, times)
    time_statistics = [volume, volume**3, μ, σ]
    time_statistics.extend(time_list)
    return time_statistics


# In[7]:


def get_statistics(f, arg, times):
    time_list = []
    for _ in range(times+1):
        t = measure_time(f, arg)
        time_list.append(t)
        
    times_without_compiling = time_list[1:]
    μ = np.mean(times_without_compiling)
    σ = np.std(times_without_compiling)
    
    return time_list, μ, σ


# In[8]:


def measure_time(f, arg):
    start_time = time.time()
    f(arg)
    finish_time = time.time()
    t = finish_time - start_time
    return t


# In[12]:


def write_to_file(f, list_of_times):
    file_name = f.__name__
    with open(f'../benchmarking/{file_name}_python.csv', 'a', newline='') as f_object:  
        writer_object = writer(f_object)
        writer_object.writerow(list_of_times)
        f_object.close()


# In[ ]:


volumes = [10,20,30,40,50,60,80,100]
times = 10
iterate_volumes(volumes, times)


# ---
