'''
This file uses the experimental data and the one obtained with the simulations using
the vertex model to show and compare the evolution of the area an perimeter of each
active cell, creating a figure: "fig_allfree.pdf", in which red (vertex model) and
orange (experiment) are area-curves, and blue (vertex model) and sky-blue (experiment)
are perimeter-curves. All curves are normalized by their values at t=0. To avoid
superposition of curves, areas and perimeters ones are displaced in the vertical axis,
similar as in the Fig.10 of the article.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


n_txt=55
n_celulas=68
dt = 0.01
d_pasos = 0.2/dt


Aexp_c = []
A_c = []
Pexp_c = []
P_c = []

texto0 = './datos_exp_celulas/0_dat_exp_cel.txt'
data0 = np.loadtxt(texto0)
A0 = data0[:,1]
P0 = data0[:,1]

def cel(c):
    Aexp_c.append([])
    A_c.append([])
    Pexp_c.append([])
    P_c.append([])
    for j in range(n_txt):
        texto1 = './datos_exp_celulas/'+str(int(j))+'_dat_exp_cel.txt'
        data1 = np.loadtxt(texto1)
        Aexp = data1[:,1]
        Pexp = data1[:,2]
        texto2 = './resultados/'+str(int(j))+'_celulas.txt'
        data2 = np.loadtxt(texto2)
        area = data2[:,3]
        per = data2[:,4]

        Aexp_c[c].append(Aexp[c])
        A_c[c].append(area[c])
        Pexp_c[c].append(Pexp[c])
        P_c[c].append(per[c])

for c in range(n_celulas):
    cel(c)
plt.close('all')
x = np.linspace(0, n_txt*(d_pasos*dt),n_txt)
pparam = dict( ylabel='$A$ ')
pparam2 = dict(xlabel='$t$', ylabel='$P$')

with plt.style.context(['science', 'grid']):
    fig, ax = plt.subplots(3, 5,figsize=(10, 6))
    fig.suptitle('Best set inner perimeter activity, $\delta=0.5$',fontsize=15)
    n = 10
    et=0

    for p in [12,14,17,19,22]: #first set of five active cells
        p1 = ax[0,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(Aexp_c[p])/Aexp_c[p][0],  '.-',markersize=5,linewidth=4, zorder=n, color='DarkOrange')
        ax[0,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(A_c[p])/A_c[p][0],  '-',markersize=2,color='red',linewidth=3,zorder=n)
        p2 = ax[0,et].plot(x, np.array(Pexp_c[p])/Pexp_c[p][0],  '.-',markersize=5,linewidth=4,zorder=n,color='dodgerblue')
        ax[0,et].plot(x, np.array(P_c[p])/P_c[p][0],  '-',markersize=2,color='mediumblue',linewidth=3,zorder=n)
        n=n-1
        ax[0,et].title.set_text('C'+str(p))
        ax[0,et].set_xlabel('$t$',fontsize=15)
        ax[0,et].xaxis.set_tick_params(labelsize=14)
        et=et+1

    et=0
    for p in [25,28,29,30,31]: #second set of five active cells
        p1 = ax[1,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(Aexp_c[p])/Aexp_c[p][0],  '.-',markersize=5,linewidth=4, zorder=n, color='DarkOrange')
        ax[1,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(A_c[p])/A_c[p][0],  '-',markersize=2,color='red',linewidth=3,zorder=n)
        p2 = ax[1,et].plot(x, np.array(Pexp_c[p])/Pexp_c[p][0],  '.-',markersize=5,linewidth=4,zorder=n,color='dodgerblue')
        ax[1,et].plot(x, np.array(P_c[p])/P_c[p][0],  '-',markersize=2,color='mediumblue',linewidth=3,zorder=n)
        n=n-1
        ax[1,et].title.set_text('C'+str(p))
        ax[1,et].set_xlabel('$t$',fontsize=15)
        ax[1,et].xaxis.set_tick_params(labelsize=14)
        et=et+1

    et=0
    for p in [32,33,37,39,44]: #third set of five active cells
        p1 = ax[2,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(Aexp_c[p])/Aexp_c[p][0],  '.-',markersize=5,linewidth=4, zorder=n, color='DarkOrange')
        ax[2,et].plot(x, 1.*np.ones(len(Aexp_c[p]))+np.array(A_c[p])/A_c[p][0],  '-',markersize=2,color='red',linewidth=3,zorder=n)
        p2 = ax[2,et].plot(x, np.array(Pexp_c[p])/Pexp_c[p][0],  '.-',markersize=5,linewidth=4,zorder=n,color='dodgerblue')
        ax[2,et].plot(x, np.array(P_c[p])/P_c[p][0],  '-',markersize=2,color='mediumblue',linewidth=3,zorder=n)
        n=n-1
        ax[2,et].title.set_text('C'+str(p))
        ax[2,et].set_xlabel('$t$',fontsize=15)
        ax[2,et].xaxis.set_tick_params(labelsize=14)

        et=et+1

    plt.setp(ax, yticks=[1,2], yticklabels=[' ', ' '],xlabel='$t$')

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('fig_allfree.pdf')
