'''
On this file we create the tissue using the experimental data from the main code,
and the classes and functions from classes_and_functions.py.
'''
#Creation of vertices
pos0 = []
vel0 = []
cells0 = []
t = []

for i in range(n_vertices):
    ady_i = []
    pos0.append((x[i],y[i],z[i]))
    vel0.append((0,0,0))
    cell_i = []
    for l in range(len(celda)):
        if v_celda[l] == i:
            cell_i.append(int(celda[l]))
    cells0.append(cell_i)
    for j in range(len(cells0[i])):
        cell = cells0[i][j]
        vs = []
        for k in range(len(celda)):
            if celda[k] == cell:
                vs.append(v_celda[k])
        if vs[0] != i and vs[len(vs)-1] != i:
            for h in range(len(vs)):
                if vs[h] == i:
                    ady_i.append([int(vs[h-1]), int(vs[h+1])])
        if vs[0] != i and vs[len(vs)-1] == i:
            ady_i.append([int(vs[len(vs)-2]), int(vs[0])])
        if vs[0] == i and vs[len(vs)-1] != i:
            ady_i.append([int(vs[len(vs)-1]), int(vs[1])])
    t.append(Vertex(tipo_v[i], pos0[i],vel0[i],cells0[i], ady_i, dt, R, KA, KP,
                     visa, visp, J))



#Creations of cells
celulas = []
for i in range(n_celulas):
    b = []
    c = []
    for j in range(len(celda)):
        if celda[j] == i:
            b.append(int(v_celda[j]))
    celulas.append(Cell(i,b,  (tc[i], delta[i], cmax[i], fracm[i]),(tcPer[i], deltaPer[i], cmaxPer[i],fracmPer[i]), dt,area00[i],per00[i], visa, visp))   #b tiene repetido b[inicial] b[final]
    b.append(b[0])


#Creation of tissue
T1 = Tissue(t, celulas, dt, R, mA,mP)
for c in borde:
    for vi in T1.cs[c].ind:
        T1.vs[vi].tipo = 0
