'''
On this file we define the classes: Vertex, Cell and Tissue, that we use in our main
code to simulate the tissue using the vertex model. We also define some useful functions
at the end.
'''
#CLASSES
class Vertex:
    def __init__(self, tipo, vert_cell, velocidad, celulas, adyacentes,
                 dt, R, KA, KP, visa, visp, J):
        self.x = vert_cell[0]
        self.tipo = tipo # 1: vertex model dinamic, 0: experimentally interpolated
        self.y = vert_cell[1]
        self.z = vert_cell[2]
        self.r = (vert_cell[0],vert_cell[1],vert_cell[2])
        self.r = (self.x, self.y, self.z)
        self.vx = velocidad[0]
        self.vy = velocidad[1]
        self.vz = velocidad[2]
        self.cells = celulas
        self.ady = adyacentes
        self.dt = dt
        self.R = R
        self.KA = KA
        self.KP = KP
        self.visa = visa
        self.visp = visp
        self.J = J
        self.x = self.x + self.vx * self.dt
        self.y = self.y + self.vy * self.dt
        self.z = self.z + self.vz * self.dt
        self.r = (self.x, self.y, self.z)
    def param(self):
        r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        return 0.5 * ( (self.R / r) - 1)
    def colapso(self):
        r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        param = self.param()
        self.x = self.x + (param * 2 * self.x)
        self.y = self.y + (param * 2 * self.y)
        self.z = self.z + (param * 2 * self.z)
        self.r = (self.x, self.y, self.z)
    def v_nuevas(self):
        visa = self.visa
        visp = self.visp
        k_p =  self.KP
        k_a =  self.KA
        J = self.J
        vx = 0
        vy = 0
        vz = 0

        for j in range(len(self.cells)):
            c = self.cells[j]
            vi = int(self.ady[j][0])
            vd = int(self.ady[j][1])

            #Proportional to perimeter variation
            l_ii = distancia(t[vi].r,self.r)
            l_id = distancia(t[vd].r,self.r)
            per0_c = per0_ec_mov[c]
            PER = 0.5*J + k_p * (celulas[c].per - per0_c)
            vx = vx +  PER * (t[vi].x - self.x)/ l_ii
            vx = vx +  PER * (t[vd].x - self.x)/ l_id
            vy = vy +  PER * (t[vi].y - self.y)/ l_ii
            vy = vy +  PER * (t[vd].y - self.y)/ l_id
            vz = vz +  PER * (t[vi].z - self.z)/ l_ii
            vz = vz +  PER * (t[vd].z - self.z)/ l_id

            #Proportional to area variation
            vec_c = celulas[c].cal_centro()
            k = div(vec_c,modulo(vec_c))
            rest = resta_vectores(t[vd].r, t[vi].r)
            v_tot = np.cross(k,rest)
            area0_c = area0_ec_mov[c]
            AREA = -0.5* k_a * (celulas[c].area - area0_c)
            vx = vx +  AREA * v_tot[0]
            vy = vy +  AREA * v_tot[1]
            vz = vz +  AREA * v_tot[2]

        return [vx, vy, vz]

class Cell(Vertex):
    def __init__(self, numero, indices, info, infoper, dt, A, P, visa, visp):
        self.ind = indices
        self.tc = info[0]
        self.delta = info[1]
        self.cmax = info[2]
        self.fracm = info[3]
        self.tcPer = infoper[0]
        self.deltaPer = infoper[1]
        self.cmaxPer = infoper[2]
        self.fracmPer = infoper[3]
        self.dt = dt
        self.n = numero
        self.area = A
        self.per = P
        self.visa = visa
        self.visp = visp

    def P_vertex(self, now):
        #medial active pulse:
        if now >= (self.tc/dt) and now <= ((self.tc+self.delta)/dt):
            valor = c_t_nuevo_gral(self.tc, self.fracm,self.cmax, self.delta, self.dt)[int(now)]
            Pulso[self.n] = valor
        #perimeter active pulse:
        if now >= (self.tcPer/dt) and now <= ((self.tcPer+self.deltaPer)/dt):
            valorPer = c_t_nuevo_gral(self.tcPer, self.fracmPer,self.cmaxPer, self.deltaPer, self.dt)[int(now)]
            PulsoPer[self.n] = valorPer

    def cal_area(self):
        v = self.ind
        area = 0
        vects = []
        for i in range(len(v)-1):
            vects.append(resta_vectores(t[int(v[i+1])].r,t[int(v[0])].r)) #HABIA ERRIR con t[i+1]
        for j in range(len(vects)-1):
            vec_c = self.cal_centro()
            k_c = div(vec_c,modulo(vec_c))
            area = area - 0.5 * np.dot(np.cross(vects[j],vects[j+1]), k_c)
        self.area = area

    def cal_per(self):
        v = self.ind
        per = 0
        aristas = []
        for i in range(len(v)-1):
            aristas.append(modulo(resta_vectores(t[int(v[i+1])].r,t[int(v[i])].r)))
        self.per = sum(aristas)

    def v_A0_nuevas(self):
        return -self.visa * (area0_ec_mov[self.n] - self.area) + (mA)*area00[self.n]+Pulso[self.n]*area00[self.n]

    def A0_nuevas(self):
        area0_ec_mov[self.n] = area0_ec_mov[self.n] + self.v_A0_nuevas() * self.dt

    def v_P0_nuevas(self):
        return -self.visp * (per0_ec_mov[self.n]-self.per)+ (mP)*per00[self.n]+PulsoPer[self.n]*per00[self.n]

    def P0_nuevas(self):
        per0_ec_mov[self.n] = per0_ec_mov[self.n] + self.v_P0_nuevas() * self.dt

    def cal_centro(self):
        c = (0,0,0)
        v = self.ind
        for i in range(len(v)-1):
            add = resta_vectores(t[int(v[i])].r,t[int(v[0])].r)
            c = suma(c,add)
        return suma(t[int(v[0])].r,div(c,len(v)-1))

class Tissue(Cell):
    def __init__(self, vertices, celulas, dt, R, mA,mP):
        self.vs = vertices
        self.cs = celulas
        self.dt = dt
        self.R = R
        self.mA = mA
        self.mP = mP
    def pulso_vertex(self, now):
        for i in range(len(self.cs)):
            self.cs[i].P_vertex(now)

    def A0_P0_nuevas_vertex(self):
        for h in range(len(self.cs)):
            self.cs[h].A0_nuevas()
            self.cs[h].P0_nuevas()

    def pos_nuevas_vertex(self,now):
        vels = []
        for i in range(len(self.vs)):
            if self.vs[i].tipo == 1:
                vels.append(self.vs[i].v_nuevas())
            else:
                vels.append([0,0,0])
        for n in range(len(self.vs)):
            if self.vs[n].tipo == 1: #positions change due to the model
                self.vs[n].vx = vels[n][0]
                self.vs[n].vy = vels[n][1]
                self.vs[n].vz = vels[n][2]

                self.vs[n].x = self.vs[n].x + self.vs[n].vx*self.vs[n].dt
                self.vs[n].y = self.vs[n].y + self.vs[n].vy*self.vs[n].dt
                self.vs[n].z = self.vs[n].z + self.vs[n].vz*self.vs[n].dt
                self.vs[n].r = (self.vs[n].x, self.vs[n].y, self.vs[n].z)

            else: #experimentally interpolated
                self.vs[n].x = evol_interpol[n][0](now)
                self.vs[n].y = evol_interpol[n][1](now)
                self.vs[n].z = evol_interpol[n][2](now)
                self.vs[n].r = (self.vs[n].x, self.vs[n].y, self.vs[n].z)

    def colaps(self):
        for m in range(len(self.vs)):
            self.vs[m].colapso()

    def cal_area_per(self):
        for i in range(len(self.cs)):
            self.cs[i].cal_area()
            self.cs[i].cal_per()

    def evol_vertex(self,now):
        self.A0_P0_nuevas_vertex()
        self.pos_nuevas_vertex(now)
        self.colaps()
        self.cal_area_per()

    def data(self):
        pos = [] #positions
        A = [] #areas
        P = [] #perimeters


        for i in range(len(self.vs)):
            pos.append((self.vs[i].x,self.vs[i].y,self.vs[i].z))
        for i in range(len(self.cs)):
            A.append(self.cs[i].area)
            P.append(self.cs[i].per)

        new01 = range(n_vertices)
        new02 = range(n_celulas)

        np.savetxt('./resultados/data_vertices.txt',  np.c_[new01, tipo_v, pos], fmt='%1.10f')
        np.savetxt('./resultados/data_celulas.txt',  np.c_[new02, area0_ec_mov,per0_ec_mov,A,P], fmt='%1.10f')



#FUNCTIONS
def distancia(vector1, vector2):
    x = vector1[0]-vector2[0]
    y = vector1[1]-vector2[1]
    z = vector1[2]-vector2[2]
    return np.sqrt(x**2 + y**2 + z**2)

def modulo(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

def resta_vectores(a,b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def suma(a,b):
    c = []
    for i in range(len(a)):
        c.append(a[i] + b[i])
    return c

def div(v,a):
    return (v[0]/a, v[1]/a, v[2]/a)

def det(a):
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]

#Active pulse function
def c_t_nuevo_gral(tcri,fracm, C1, delta, dt):
    tm = tcri+delta*fracm
    delta1 = (tm-tcri)
    delta2 = (delta+tcri-tm)
    time1 = np.arange(0,delta1, dt)
    time2 = np.arange(0,delta2, dt)
    c = []
    for i in range(int(len(time1))):
        c.append(-C1*np.sin((np.pi/delta1)*(time1[i])))
    for i in range(int(len(time2))):
        c.append(C1*(1*delta1/delta2)*np.sin((np.pi/delta2)*(time2[i])))
    pasos_1hr = 1./(d_pasos*dt)
    c1 = np.zeros(int(pasos_1hr * tcri * d_pasos))
    total_pasos = total_txt*d_pasos
    c2 = np.zeros(int(total_pasos-len(c)-len(c1)))
    c3 = np.concatenate((np.concatenate((c1,c), axis=0),c2), axis=0)
    return c3

def calcular_E_area():
    E_tot = 0
    for i in range(total_txt):
        E_i = 0
        area_sim = np.loadtxt('./resultados/'+str(i)+'_celulas.txt')[:,3]
        area_exp = np.loadtxt('./datos_exp_celulas/' +str((i))+'_dat_exp_cel.txt')[:,1]
        for j in range(n_celulas):
            delta = 0
            for number in borde:
                if number == j:
                    delta = delta + 1
            if delta == 0:
                E_i = E_i + (area_sim[j] - area_exp[j])**2
        E_tot = E_tot + E_i
    return E_tot /(total_txt*(n_celulas-len(borde))*(Aprom**2))

def calcular_E_per():
    E_tot = 0
    for i in range(total_txt):
        E_i = 0
        P_sim = np.loadtxt('./resultados/'+str(i)+'_celulas.txt')[:,4]
        P_exp = np.loadtxt('./datos_exp_celulas/' +str((i))+'_dat_exp_cel.txt')[:,2]
        for j in range(n_celulas):
            delta = 0
            for number in borde:
                if number == j:
                    delta = delta + 1
            if delta == 0:
                E_i = E_i + (P_sim[j] - P_exp[j])**2
        E_tot = E_tot + E_i
    return E_tot/(total_txt*(n_celulas-len(borde))*(Pprom**2))
