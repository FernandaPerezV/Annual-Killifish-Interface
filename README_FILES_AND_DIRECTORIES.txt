#####  MainFile.py  #####
This code produces the evolution of the epithelial tissue when using the best
parameters (global and for the active pulses). On this file we define the numerical
parameters of the simulation. Then, we import the experimental data. The program runs "classes_and_functions.py",
Which has all the definitions of classes and functions. Later, we create the tissue (the program runs other file to do that:
tissue_creation.py). Finally, the dynamic starts, and we save information
every d_pasos steps, about the vertices and cells, on the directory "resultados".
By the end, we calculate the optimisation functions related to the areas and perimeters,
E1 and E2, that use the experimental information from the directory "datos_exp_celulas". Finally, E1 and E2 and are printed. 


#####  classes_and_functions.py  #####
On this file we define the classes: Vertex, Cell and Tissue, that we use in our main
code to simulate the tissue using the vertex model. We also define some useful functions
at the end.

#####  tissue_creation.py  #####
On this file we create the tissue using the experimental data from the main code,
and the classes and functions from pre_ep1_interpol.py.


#####  resultados  #####
In this directory we save information from the simulations about the vertices and cells.


#####  evlVecinosDeVertice_m0  #####
col0: ID of a vertex
col1: ID of a vertex that is joined to the one of the col0 by an cellular edge.
The first ID is 1.


#####  celdas_mauricio  #####
col0: ID of a cell
col1: ID of a vertex that belongs to the cell of the col0. These ID's are
clock-wise ordered.
The first ID is 1.


#####  pos_time_sindrift_esfera.txt  #####
col0: ID of a vertex
col1: experimental time step
col2: x position of the vertex of the col0
col2: y position of the vertex of the col0
col2: -z position of the vertex of the col0
The first ID is 1.


#####  pos0  #####
First 134 lines of the pos_time_sindrift_esfera.txt file, i.e., experimental
information at the beginning of the experiment.
The first ID is 1.


#####  APinicial_m_3.txt  #####
The first ID is 0.
col0: ID of a cell
col1: area of the cell of the col0, calculated by the triangularization method,
using the positions of pos_time_sindrift_esfera.txt, and the topology given by
celdas_mauricio.
col2: perimeter of the cell of the col0, calculated using the euclidian distances
between the positions of pos_time_sindrift_esfera.txt, and the topology given by
celdas_mauricio.


#####  pulso.txt  #####
col0: ID of an medial-active cell
col1: Time of the start of the active pulse.
col2: Duration of the active pulse.
col3: 1 if the pulse exists, 0 if the pulse does not exist.
col4: C1, amplitude of the destruction phase of the active pulse.
col5: Fraction between the destruction and construction phases of the active pulse.


#####  pulsoPer.txt  #####
Same as in pulso.txt, but now for perimeter activity.


#####  datos_exp_celulas  #####
This directory contains one file per experimental step, expstep_data_exp_cel.txt,
in which:
col0: ID of a cell
col1: area of the cell of the col0, calculated by the triangularization method,
using the positions of pos_time_sindrift_esfera.txt, and the topology given by
celdas_mauricio.
col2: perimeter of the cell of the col0, calculated using the euclidian distances
between the positions of pos_time_sindrift_esfera.txt, and the topology given by
celdas_mauricio.
col3: x position of the center of the cell of the col0.
col4: y position of the center of the cell of the col0.
The first ID is 0.


#####  plot_AP_final.py  #####
This file uses the experimental data and the one obtained with the simulations using
the vertex model to show and compare the evolution of the area an perimeter of each
active cell, creating a figure: "fig_allfree.pdf", in which red (vertex model) and
orange (experiment) are area-curves, and blue (vertex model) and sky-blue (experiment)
are perimeter-curves. All curves are normalized by their values at t=0. To avoid
superposition of curves, areas and perimeters ones are displaced in the vertical axis,
similar as in the Fig.10 of the article.

#####  fig_allfree.pdf  #####
Figure created when running plot_AP_final.py, described previously.