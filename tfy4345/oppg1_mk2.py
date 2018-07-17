# oppgave 2, klassisk mekanikk
# Stian Hartman, Erik Liodden

import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import time


class Planet:
    def __init__(self, name, mass, r0, v0):
        self.name = name
        self.m = mass
        self.r = np.array(r0)
        self.v = np.array(v0)
        self.X = np.zeros(N)
        self.Y = np.zeros(N)
        self.number_of_updates = 0


def update_v2(planet):
    F = force_v2(planet)
    # v_halv = np.zeros((np.size(planet),2))
    for i in range(np.size(planet)):
        v_halv = planet[i].v + F[i] * planet[i].r * (dt / 2)
        planet[i].r = planet[i].r + v_halv * dt
        F = force_v2(planet)
        planet[i].v = v_halv + F[i] * planet[i].r * (dt / 2)
        planet[i].X[planet[i].number_of_updates] = planet[i].r[0]
        planet[i].Y[planet[i].number_of_updates] = planet[i].r[1]
        planet[i].number_of_updates += 1
    return planet

def potential_v2(planet):
    V = 0
    for i in range(np.size(planet)):
        V += - G * planet[i].m * M / np.linalg.norm(planet[i].r)
    for i in range(np.size(planet)-1):
        j = i + 1
        while j < np.size(planet):
            V += - G * planet[i].m * planet[j].m / np.linalg.norm(planet[i].r - planet[j].r)
            j += 1
    return V

def force_v2(planet):
    # scalar version of force
    F = np.zeros(np.size(planet))
    for i in range(np.size(planet)):
        F[i] = - G * M / (np.linalg.norm(planet[i].r) ** 3)
    for i in range(np.size(planet)):
        j = 0
        while j < np.size(planet):
            if i != j:
                F[i] += - G * planet[i].m * planet[j].m / np.linalg.norm(planet[i].r - planet[j].r)**3
            j += 1
    return F

def total_energy(planet):
    pot = potential_v2(planet)
    E_tot = 0
    for i in range(np.size(planet)):
        E_kin = 1/2 * planet[i].m * np.linalg.norm(planet[i].v)**2
        E_tot += E_kin
    return E_tot + pot

def plotter_1(planet):
    for i in range(np.size(planet)):
        plt.plot(planet[i].X[:planet[i].number_of_updates], planet[i].Y[:planet[i].number_of_updates])
        plt.scatter(0, 0, s = 30)
        plt.scatter(planet[i].X[planet[i].number_of_updates-1], planet[i].Y[planet[i].number_of_updates-1], s = 10*planet[i].m)



# diverse globale variabler
dt = 0.001
N = 2000
G = 1  # gravitational constant
M = 750  # solar mass
energy = np.zeros(N)

# lager planeter og setter initialverdiene

nr1 = Planet('nr1', 1, [4, 0], [0, 10])
nr2 = Planet('nr2', 2, [3, 0], [0, 13])
# nr3 = Planet('nr3', 3, [0, 5], [-10, 0])
# nr4 = Planet('nr4', 4, [5, 1], [0, -13])
# nr5 = Planet('nr5', 4, [-2, -2], [10, -2])
# nr6 = Planet('nr6', 8, [0, -8], [10, 0])

# planets = np.array([nr1, nr2, nr3, nr4, nr5, nr6])
planets = np.array([nr1, nr2])
plt.ion()


for i in range(N):
    planets = update_v2(planets)
    energy[i] = total_energy(planets)

    if i % 10 == 0:
        plt.cla()
        plt.xlim([-10, 10])
        plt.ylim([-10, 10])
        plotter_1(planets)
        plt.draw()
        time.sleep(0.005)
    if i % (N/100) == 0:
        print('%i%% complete' % (i/N*100))


plt.figure()
plt.ioff()
plt.plot(np.arange(0, N*dt, dt), energy)
plt.show()
