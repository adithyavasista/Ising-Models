import numpy as np
import random
import time
import math

p = 500 # no. of time slices which is the Trotter number y axis
n = 500 # no. of spins on the horizontal axis
J = 1
Jx = J / p  # interaction energy on horizontal axis
temp = 0.06 # the temperature
min_h = 0.02
max_h = 6
h_int = 0.02
h = min_h
accept = 0

def initialize(s):
    for i in range(p):
        for j in range(n):
            s[i][j] = 1  # uniform start

def metropolis(s, new_energy, Jy):
    global accept
    k = random.randint(0, n - 1)  # choose site k on x axis for flipping
    l = random.randint(0, p - 1)  # choose site l on y axis for flipping
    r = random.random()
    delta_E = 2 * s[l][k] * (Jx * (s[l][(k - 1) % n] + s[l][(k + 1) % n]) + Jy * (s[(l - 1) % p][k] + s[(l + 1) % p][k]))
    if math.exp(-delta_E / temp) > r:
        s[l][k] *= -1
        new_energy += delta_E
        accept += 1
        return True
    else:
        return False

def one_mc(s, new_energy, Jy):
    for i in range(p * n):
        if metropolis(s, new_energy, Jy):
            global accept
            accept += 1

def mag_per_spin(s):
    sum_spins = np.sum(s)
    return sum_spins / (n * p)

def energy(s, Jy):
    sum_energy = 0.0
    for i in range(p):
        for j in range(n):
            sum_energy += Jx * s[i][j] * s[i][(j + 1) % n] + Jy * s[i][j] * s[(i + 1) % p][j]
    term = (0.5 * math.sinh(2 * h / (temp * p))) ** ((p * n) / 2.0)
    if term > 0:
     sum_energy += math.log(term)
    else:
     sum_energy += 0  # or some other value that makes sense in your context
    return -1 * sum_energy / (n * p)

with open("2D_Ising_Trot_mag.dat", "w") as outfile:
    outfile.write("field\tmag_per_spin\n")
    mcsteps = 120  # no. of Monte Carlo steps
    mag_sum = 0
    old_energy = 0

    s = np.ones((p, n), dtype=int)  # initialize the spin array
    initialize(s)
    Jy = 0.5 * temp * math.log(1 / math.tanh(min_h / (temp * p)))
    old_energy = energy(s, Jy)
    new_energy = old_energy

    while h < max_h:
        Jy = 0.5 * temp * math.log(1 / math.tanh(h / (temp * p)))
        accept = 0
        random.seed(time.time())
        therm_steps = int(0.2 * mcsteps)
        for _ in range(therm_steps):
            one_mc(s, new_energy, Jy)

        mag_sum = 0.0
        for _ in range(mcsteps):
            one_mc(s, new_energy, Jy)
            mag_sum += mag_per_spin(s)

        outfile.write(f"{h:.2f}\t{mag_sum / mcsteps:.6f}\n")
        h += h_int
