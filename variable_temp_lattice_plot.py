import numpy as np
import matplotlib.pyplot as plt
# Parameters
L = int(input("Enter the number of lattice points in one dimension: "))
niter = int(input("Enter the number of iterations: "))
T = 2.0
J_ising = 1.0
M=0
# Initialize lattice with random spins
spin = np.zeros((L,L),dtype=int)
for i in range(L):
    for j in range(L):
        a=np.random.rand()
        if a<0.5:
            spin[i,j]=1
        else:
            spin[i,j]=-1
        M=M+ spin[i,j]
storage=np.zeros(niter)
counter=np.zeros(niter)
# Functions to compute neighbors with periodic boundary conditions (PBC)
def get_neighbors(i, j, L):
    a = (i + 1) 
    b = (i - 1) 
    c = (j + 1) 
    d = (j - 1) 
    if i==(L-1):
        a=0
    if i==0:
        b=L-1
    if j==0:
        d=L-1
    if j==(L-1):
        c=0
    return a, b, c, d
# Initialize energy and magnetization
E = 0.0
N = L * L

for i in range(L):
    for j in range(L):
        a, b, c, d = get_neighbors(i, j, L)
        E -= J_ising * spin[i, j] * (spin[a, j] + spin[b, j] + spin[i, c] + spin[i, d])

mag = M / float(N)
mavg=0

print(f'Initial energy E, E per spin = {E}, {E / float(N)}')
print(f'Initial magnetization M, M per spin = {M}, {mag}')
# Spin flipping and Metropolis algorithm
t=15
for t in range(10,20):
    T=t/10.0
    for time in range(niter):
        for p in range(L):
            for k in range(L):
                E_1 =0
                E_2=0
                i = np.random.randint(L)
                j = np.random.randint(L)
                a, b, c, d = get_neighbors(i, j, L)
                E_2 -= J_ising * spin[i, j] * (spin[a, j] + spin[b, j] + spin[i, c] + spin[i, d])
                spin[i, j] *= -1
                E_1 -= J_ising * spin[i, j] * (spin[a, j] + spin[b, j] + spin[i, c] + spin[i, d])
                dE =E_1-E_2
                if dE <= 0: 
                    E = dE+ E
                    M += 2 * spin[i, j]
                if dE >0:
                    r=np.random.rand()
                    if r < np.exp(-dE / T):
                         E += dE
                         M += 2 * spin[i, j]
                    else :
                        spin[i, j] *= -1
        if(time<1000):       
            mavg=0
        else:
            mavg=abs(M)+mavg


        mag = M/ float(N)
        eng=E/float(N)
        storage[time]=mag
        counter[time]=time
    print(f'Final energy E, E per spin = {E}, {E / float(N)}')
    print(f'Final magnetization M, M per spin = {M}, {mag}')
    mavg=mavg/(float(N)*float(niter-1000))
    print(f'average magnetization  M per spin = {mavg}')
    plt.plot(range(niter), storage, label=f'T={T}')
plt.xlabel("Iterations")
plt.ylabel("Magnetization per spin")
plt.title("Magnetization per spin for various temperatures")
plt.legend()
plt.show()