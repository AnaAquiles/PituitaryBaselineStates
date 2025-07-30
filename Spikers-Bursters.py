
import numpy as np
import matplotlib.pyplot as plt
import random

"""
         The following functions was developped by Mehran Fazli and Richard Bertram
         available here: https://www.math.fsu.edu/~bertram/software/pituitary/ 

         and published in the following DOI:10.1371/journal.pcbi.1011811 

         Disclaimer: The following version contain adaptation and editions published in: 
        Results of this modification are showed in Figure 2 from the same review.

"""
plt.style.use('fivethirtyeight')

def ode_solver(IC,Nt,ggj,gbk,mat):

    vn=-5
    kc=0.12
    ff=0.005
    vca=60
    vk=-75
    vl = -50.0
    gk=2.5
    cm=5
#    gbk=1
    gca=2.1
    gsk=2
    vm=-20
    vb=-5
    sn=10
    sm=12
    sbk=2
    taun=30
    taubk=5
    ks=0.4
    alpha=0.0015
    gl=0.2
#    ggj=0.05
    dt=0.5

    #//////////////////////////////////
    # LOOP OF INITIAL CONDITION



    vtim=np.zeros((num_cell,Nt))
    ntim=np.zeros((num_cell,Nt))
    ctim=np.zeros((num_cell,Nt))
    btim=np.zeros((num_cell,Nt))

    vv=np.zeros(num_cell)
    nn=np.zeros(num_cell)
    cc=np.zeros(num_cell)
    bb=np.zeros(num_cell)
    

    for i in range(0,num_cell):
        vv[i]= IC[0,i]
        nn[i]= IC[1,i]
        cc[i]= IC[2,i]
        bb[i]= IC[3,i]


    for t in range(0,Nt): #/////////////////////////////////////////////////SOLVER TIME-STEP ITERATION

        ninf=0
        bkinf=0
        minf=0
        cinf=0
        ica=0
        isk=0
        ibk=0
        ikdr=0
        il=0

        k1v=0
        k2v=0
        k3v=0
        k4v=0
        k1n=0
        k2n=0
        k3n=0
        k4n=0
        k1c=0
        k2c=0
        k3c=0
        k4c=0
        k1b=0
        k2b=0
        k3b=0
        k4b=0
        #time.sleep(0.0001)
    
        vv_old=np.zeros(num_cell)
        for i in range(0,num_cell):
            vv_old[i]=vv[i]


        for i in range(0,num_cell): #/////////////////////// LOOP IN CELLS

            #//////////////////////////////////////////////////////////////////////////// FIRST STEP RK
            ninf=1/(1+np.exp((vn-vv[i])/sn))
            bkinf=1/(1+np.exp((vb-vv[i])/sbk))
            minf=1/(1+np.exp((vm-vv[i])/sm))
            cinf=((cc[i])**2)/(((cc[i])**2)+ks*ks)

            ica=gca*minf*(vv[i]-vca)
            isk=gsk*cinf*(vv[i]-vk)
            ibk=gbk[i]*bb[i]*(vv[i]-vk)
            ikdr=gk*nn[i]*(vv[i]-vk)
            il = gl*(vv[i]-vl)

            igj=0;
            for h in range(0,num_cell):
                igj=igj+mat[i][h]*ggj*(vv[i]-vv_old[h])

            k1v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k1n = dt*((ninf-nn[i])/taun)
            k1c = dt*(-ff*(alpha*ica+kc*cc[i]))
            k1b = dt*(-(bb[i]-bkinf)/taubk)
            #//////////////////////////////////////////////////////////////////////////// SECOND STEP RK
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k1v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k1v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k1v))/sm))
            cinf=((cc[i] + 0.5*k1c)**2)/(((cc[i] + 0.5*k1c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k1v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k1v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k1b)*((vv[i] + 0.5*k1v)-vk)
            ikdr=gk*(nn[i] + 0.5*k1n)*((vv[i] + 0.5*k1v)-vk)
            il = gl*((vv[i] + 0.5*k1v)-vl)
            igj=0
            for h in range(0,num_cell):
                igj=igj+mat[i][h]*ggj*((vv[i] + 0.5*k1v)-vv_old[h])

            k2v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k2n = dt*((ninf-(nn[i] + 0.5*k1n))/taun)
            k2c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k1c)))
            k2b = dt*(-((bb[i] + 0.5*k1b)-bkinf)/taubk)
            #//////////////////////////////////////////////////////////////////////////////// THIRD STEP RK
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k2v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k2v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k2v))/sm))
            cinf=((cc[i] + 0.5*k2c)**2)/(((cc[i] + 0.5*k2c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k2v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k2v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k2b)*((vv[i] + 0.5*k2v)-vk)
            ikdr=gk*(nn[i] + 0.5*k2n)*((vv[i] + 0.5*k2v)-vk)
            il = gl*((vv[i] + 0.5*k2v)-vl);
            igj=0;
            for h in range(0,num_cell):
                igj=igj+mat[i][h]*ggj*((vv[i] + 0.5*k2v)-vv_old[h])

            k3v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k3n = dt*((ninf-(nn[i] + 0.5*k2n))/taun)
            k3c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k2c)))
            k3b = dt*(-((bb[i] + 0.5*k2b)-bkinf)/taubk)
            #////////////////////////////////////////////////////////////////////////////////// FOURTH STEP RK
            ninf=1/(1+np.exp((vn-(vv[i] + 0.5*k3v))/sn))
            bkinf=1/(1+np.exp((vb-(vv[i] + 0.5*k3v))/sbk))
            minf=1/(1+np.exp((vm-(vv[i] + 0.5*k3v))/sm))
            cinf=((cc[i] + 0.5*k3c)**2)/(((cc[i] + 0.5*k3c)**2)+ks*ks)

            ica=gca*minf*((vv[i] + 0.5*k3v)-vca)
            isk=gsk*cinf*((vv[i] + 0.5*k3v)-vk)
            ibk=gbk[i]*(bb[i] + 0.5*k3b)*((vv[i] + 0.5*k3v)-vk)
            ikdr=gk*(nn[i] + 0.5*k3n)*((vv[i] + 0.5*k3v)-vk)
            il = gl*((vv[i] + 0.5*k3v)-vl)
            igj=0;
            for h in range(0,num_cell):
                igj=igj+mat[i][h]*ggj*((vv[i] + 0.5*k3v)-vv_old[h])

            k4v = dt*(-(ica+isk+ibk+ikdr+igj+il)/cm)
            k4n = dt*((ninf-(nn[i] + 0.5*k3n))/taun)
            k4c = dt*(-ff*(alpha*ica+kc*(cc[i] + 0.5*k3c)))
            k4b = dt*(-((bb[i] + 0.5*k3b)-bkinf)/taubk)
            #////////////////////////////////////////////////////////////////////////////////// FINAL STEP RK

            vv[i] = vv[i] + (1.0/6.0)*(k1v + 2*k2v + 2*k3v + k4v)
            nn[i] = nn[i] + (1.0/6.0)*(k1n + 2*k2n + 2*k3n + k4n)
            cc[i] = cc[i] + (1.0/6.0)*(k1c + 2*k2c + 2*k3c + k4c)
            bb[i] = bb[i] + (1.0/6.0)*(k1b + 2*k2b + 2*k3b + k4b)
            vtim[i,t]=vv[i]
            ntim[i,t]=nn[i]
            ctim[i,t]=cc[i]
            btim[i,t]=bb[i]

    return vtim, ntim, ctim, btim

"""
      Create a new binary interaction with choosing the level of population of spikers
      and bursters with changing the variable [probability_of_one], default value = 0.8

"""
def random_adjacency_matrix(n, probability_of_one=0.8):
    matrix = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):  # only fill upper triangle to avoid redundancy
            value = 1 if random.random() < probability_of_one else 0
            matrix[i][j] = value
            matrix[j][i] = value  # ensure symmetry

    return matrix



ToyMat2 = np.array(random_adjacency_matrix(10))
RavelMat2 = ToyMat2.ravel()

ggj=0.05        # medium coupling strength
gbk = ([1,0])

dt=0.5
sec=60
Ntim=int((sec*1000*(1/dt)))
num_cell=len(ToyMat2)   ### matrix
IC=np.zeros((4,num_cell))

for i in range(0,num_cell):
    IC[0,i]= -60
    IC[1,i]= 0.1
    IC[2,i]= 0.1
    IC[3,i]= 0.1

####          Perform and initializate the simulation
###       gbk - is the array of bursters #1 and spiker #0
####   Be sure to get your ggj, gbk and SN_sim(adjacency matrix) ready

vtraj, ntraj, ctraj, btraj=ode_solver(IC,Ntim,ggj,RavelMat2, ToyMat2)

Ravelmat = RavelMat2

"""
     Visualize the results, usign the same function plots from 
                 Mehran Fazli and Richard Bertram

                 Revisited version 
"""

color_opt=['tab:red','tab:blue']
cell_type_opt=['burster', 'spiker']

gbk_celltype=[1,0]
fig=plt.figure(figsize=(18,12))

cols=list(islice(cycle(color_opt),len(Ravelmat)))
cell_type=list(islice(cycle(cell_type_opt),len(Ravelmat)))


ax = fig.add_subplot(3,1,1)
xx=np.arange(0,sec, 0.0005)
for i in range(num_cell):
    plt.plot(xx,vtraj[i,0:Ntim], c=cols[i], lw=4)
plt.xlabel('time (sec)')
plt.ylim([-70, 20])
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)



ax = fig.add_subplot(3,1,2)
for i in range(num_cell):
   
    plt.plot(xx,ctraj[i,0:Ntim], c=cols[i], lw=4)
plt.ylim([0.25, 0.37])
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)

ax = fig.add_subplot(3,1,3)
alp=5
half=0.6
for i in range(num_cell):
    plt.plot(xx,1/(1+np.exp(-alp*((ctraj[i,0:Ntim]-0.27)/(0.352-0.27)-half))), c=cols[i], lw=4)

plt.xlabel('time (sec)', fontsize=25)
plt.ylim([0, 1])
plt.yticks([0.2, 0.5, 0.8])

ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.setp(ax.get_xticklabels(), fontsize=25)
plt.setp(ax.get_yticklabels(), fontsize=25)


plt.tight_layout()
plt.show()


fig, ax = plt.subplots()
im = ax.imshow(vtraj,aspect='auto')
ax.set_title('Population Vm')


plt.show()

fig, ax = plt.subplots()
im = ax.imshow(ctraj,aspect='auto')
ax.set_title('Population Calcium')


plt.tight_layout()
plt.show()

xx=np.arange(0,sec, 0.0005)
for i in range(num_cell):
    fig=plt.figure(figsize=(18,12))

    ax = fig.add_subplot(1,1,1)


    plt.plot(xx,vtraj[i,0:Ntim], c=cols[i], lw=4)
    plt.xlabel('time (sec)')
    plt.ylim([-70, 20])
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.tight_layout()
    plt.show()
    plt.close(fig)



alp=5
half=0.6
for i in range(num_cell):
    fig=plt.figure(figsize=(18,12))
    ax = fig.add_subplot(1,1,1)
    plt.plot(xx,1/(1+np.exp(-alp*((ctraj[i,0:Ntim]-0.27)/(0.352-0.27)-half))), c=cols[i], lw=4)

    plt.xlabel('time (sec)', fontsize=25)
    plt.ylim([0, 1])
    plt.yticks([0.2, 0.5, 0.8])
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
plt.tight_layout()
plt.show()
