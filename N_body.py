import numpy as np
import scipy.linalg as la
from scipy.linalg import sqrtm

def pauli_matrices():
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_y = np.array([[0, -1j], [1j, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])
    return sigma_x, sigma_y, sigma_z

def create_ising_hamiltonian(J, h, N):
    sigma_x, _, sigma_z = pauli_matrices()
    I = np.eye(2)
    
    H = np.zeros((2**N, 2**N), dtype=complex)
    
    # Transverse field term
    for i in range(N):
        Hz = h * np.kron(np.kron(np.eye(2**i), sigma_z), np.eye(2**(N-i-1)))
        H += Hz
        
    # Interaction term
    for i in range(N-1):
        Hx = -J * np.kron(np.kron(np.kron(np.eye(2**i), sigma_x), sigma_x), np.eye(2**(N-i-2)))
        H += Hx
    
    return H
# Example usage
N = 13# Total number of qubits
J = 1.0  # Interaction strength
# for o in range(1,200):
for o in range(1,400):
    h1=o/100.0
    sigma_y = np.array([[0, -1j], [1j, 0]])
    H = create_ising_hamiltonian(J, h1, N)
    #print(H)
    eigenvalues, eigenvectors = la.eigh(H)
    #print(eigenvectors)
    #print(eigenvalues)
    eigenvalues, eigenvectors = la.eigh(H)
    state_vector = eigenvectors[:, 0]
    state=state_vector.T.reshape(-1, 1)
    rho=(state_vector*state)
    rho_reduced = np.zeros((4, 4), dtype=complex)
    c = 0
    counter = 0
    counterj = 2 ** (N - 2)
    rho_reduced = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        r = counterj
        for j in range(4):
            if (i <= j):
                if (i == j):
                    for w in range(2 ** (N - 2)):
                        rho_reduced[j, i] = rho_reduced[j, i] + rho[c, c]
                        
                        c = c + 1
                    rho_reduced[i, j] = rho_reduced[j, i]
                else:
                    col = counter
                    for w in range(2 ** (N - 2)):
                        rho_reduced[j, i] = rho_reduced[j, i] + rho[r, col]
                    
                        col = col + 1
                        r = r + 1
                    rho_reduced[i, j] = rho_reduced[j, i]
        counter = counter + 2 ** (N - 2)
        counterj = counterj + 2 ** (N - 2)
    rho_tilda=np.kron(sigma_y, sigma_y)@(rho_reduced@(np.kron(sigma_y, sigma_y)))
    R=sqrtm(sqrtm(rho_reduced)@rho_tilda@sqrtm(rho_reduced))
    R=np.real(R)
    eigenvalues = np.linalg.eigvals(R)
    eigenvalues=np.sort(eigenvalues)
    C=np.max([0,(eigenvalues[3]-eigenvalues[2]-eigenvalues[1]-eigenvalues[0])])
    print(C)           
            
            

            
        




























































































            
    




                    

















