'''
This code produces the evolution of the epithelial tissue when using the best
parameters (global and for the active pulses). On this file we define the
numerical parameters of the simulation. Then, we import the experimental data.
The program runs "classes_and_functions.py", which has all the definitions of
classes and functions. Later, we create the tissue (the program runs other file
to do that: tissue_creation.py). Finally, the dynamic starts, and we save
information every d_pasos steps, about the vertices and cells, on the directory
"resultados". By the end, we calculate the optimisation functions related to the
areas and perimeters, E1 and E2, that use the experimental information from the
directory "datos_exp_celulas". Finally, E1 and E2 and are printed.
'''
from __future__ import division
import numpy as np
from scipy import interpolate
from shutil import copy

#####----------------- EXPERIMENTAL AND MODEL PARAMETERS -----------------#####
[visp, visa, KP, KA, mP, mA,J] = [0., 3.2,4.82, 0.0025, 0.0, 0.051,0]
#model parameters
dt = 0.01 #simulation's integration time step
R = 593.566615 #experimental egg's radious
n_celulas=68 #total amount of cells
n_vertices = 134 #total amount of vertices
experimental_dt = 0.2 #the experimental data are taken every 0.2h
d_pasos = experimental_dt/dt #amount of simulation's time steps needed to reach
#experimental_dt
total_txt = 55 #total amount of times in which data is saved (total amount of
#.txt files for vertices and cells
t_mucho_mas_final = total_txt + 5 #any time larger than total_txt. Is the preset
#starting time for the active pulse of every cell of the tissue
borde = [26,35,45,52,58,64,65,67,66,61,59,55,46,34,27,18,9,5,0,2,1,3,4,6,10,
         16,60,62,63] #cells of the border
Aprom = 7658.7728 #mean cellular area in t=0
Pprom = 334.5 #mean cellular perimeter in t=0
#####---------------------------------------------------------------------#####

#####------------------------ EXPERIMENTAL DATA --------------------------#####
#Experimental data of vertices and cells
datos_m = np.loadtxt('evlVecinosDeVertice_m0')-1
vert_m = datos_m[:,0]
vecino_m = datos_m[:,1]
celdas = np.loadtxt('celdas_mauricio')-1
celda = celdas[:,0]
v_celda = celdas[:,1]
data = np.loadtxt('pos0')
x = data[:,2]
y = data[:,3]
z = -1*data[:,4]
tipo_v = data[:,1]

#Interpolated data
data_inter = np.loadtxt('pos_time_sindrift_esfera.txt')
datax = data_inter[:,2]
datay = data_inter[:,3]
dataz = -1*data_inter[:,4]
evol_interpol = []
tiempo_interpol = np.linspace(0,d_pasos*total_txt*dt, total_txt+1)
for i in range(n_vertices):
    dx = []
    dy = []
    dz = []
    for j in range(total_txt+1):
        pos_vert_time = i+n_vertices*j
        dx.append(datax[pos_vert_time])
        dy.append(datay[pos_vert_time])
        dz.append(dataz[pos_vert_time])
    fx = interpolate.interp1d(tiempo_interpol, dx)
    fy = interpolate.interp1d(tiempo_interpol, dy)
    fz = interpolate.interp1d(tiempo_interpol, dz)
    evol_interpol.append([fx,fy,fz])

#creation of arrays for the simulation
APinicial = np.loadtxt('APinicial_m_3.txt')
area0_vis = np.copy(APinicial[:,1])
area0_ec_mov = np.copy(APinicial[:,1])
area00 = np.copy(APinicial[:,1])
per0_vis = np.copy(APinicial[:,2])
per0_ec_mov = np.copy(APinicial[:,2])
per00 = np.copy(APinicial[:,2])
#####---------------------------------------------------------------------#####

#####--------------------- MEDIAL PULSES INFORMATION ---------------------#####
#initial assumpsions
Pulso = np.zeros(n_celulas)
tc = t_mucho_mas_final * np.ones(n_celulas) #initially, all the cells have a
#starting time for and active pulse at t_mucho_mas_final
delta =  np.ones(n_celulas) #initially, all the cells have an active pulse that
#last 1
cmax =  np.ones(n_celulas) #initially, all the cells have an active pulse with
#a destruction amplitud of 1
fracm = 0.5*np.ones(n_celulas) #initially, all the cells have an active pulse
#with same time for destruction and for creation

#information of the pulses we will actually see in the simulation
datos_pulsos = np.loadtxt('pulso.txt')
for i in range(len(datos_pulsos[:,0])):
    celdaarea = datos_pulsos[:,0][i] #number (id) of the active cell
    tc[int(celdaarea)] = datos_pulsos[:,1][i] #active pulse's starting time
    delta[int(celdaarea)] = datos_pulsos[:,2][i] #total active time
    cmax[int(celdaarea)] = datos_pulsos[:,3][i] * datos_pulsos[:,4][i]
    #destruction amplitude
    fracm[int(celdaarea)] = datos_pulsos[:,5][i] #destruction time/creation time
#####---------------------------------------------------------------------#####

#####-------------------- PERIMETER PULSES INFORMATION -------------------#####
PulsoPer = np.zeros(n_celulas)
tcPer = t_mucho_mas_final * np.ones(n_celulas)
deltaPer = np.ones(n_celulas)
cmaxPer = np.ones(n_celulas)
fracmPer = 0.5*np.ones(n_celulas)
datos_pulsosPer = np.loadtxt('pulsoPer.txt')

for i in range(len(datos_pulsosPer[:,0])):
    celdaPer = datos_pulsosPer[:,0][i]
    tcPer[int(celdaPer)] = datos_pulsosPer[:,1][i]
    deltaPer[int(celdaPer)] = datos_pulsosPer[:,2][i]
    cmaxPer[int(celdaPer)] = datos_pulsosPer[:,3][i] * datos_pulsosPer[:,4][i]
    fracmPer[int(celdaPer)] = datos_pulsosPer[:,5][i]
#####---------------------------------------------------------------------#####

#Functions and classes: vertices-cells-tissue
exec(open("classes_and_functions.py").read())

#Creation of the tissue
exec(open("tissue_creation.py").read())

#Main program:
T1.data()
T1.cal_area_per()
copy('./resultados/data_vertices.txt','./resultados/0_vertices.txt')
copy('./resultados/data_celulas.txt','./resultados/0_celulas.txt')

for i in range(int(d_pasos * total_txt)):
    T1.pulso_vertex((i+1))
    T1.evol_vertex((i+1)*dt)
    c = (i + 1) % d_pasos
    d = int((i + 1) / d_pasos)
    if c == 0 :
        T1.data()
        copy('./resultados/data_vertices.txt', './resultados/'+ str(d) + '_vertices.txt')
        copy('./resultados/data_celulas.txt', './resultados/'+str(d) + '_celulas.txt')

E1 = calcular_E_area()
E2 = calcular_E_per()
print(E1, E2)
