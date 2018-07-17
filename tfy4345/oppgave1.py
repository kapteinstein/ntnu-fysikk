# oppgave 1, klassisk mekanikk
# Stian Hartman, Erik Liodden

import numpy as np
import matplotlib. pyplot as plt


def plotter():
    plt.clf()
    plt.subplot(121)
    plt.plot(X, Y)
    plt.scatter(0, 0, s = 30)
    plt.xlim([-0.5, 1.5])
    plt.ylim([-1, 1])
    plt.title('Numerically calculated orbit')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.subplot(122)
    plt.plot(X_eksakt, Y_eksakt)
    plt.scatter(0, 0, s = 30)
    plt.xlim([-0.5, 1.5])
    plt.ylim([-1, 1])
    plt.title('Exact orbit')
    plt.xlabel("x")
    plt.ylabel("y")

    plt.figure()
    ax = plt.subplot(121)
    ax.plot(Tid, E_kin/np.abs(E_kin[0]+E_pot[0]), "b", label = "E_kinetic")
    ax.plot(Tid, E_pot/np.abs(E_kin[0]+E_pot[0]), "r", label = "E_potential")
    ax.plot(Tid, (E_kin+E_pot)/np.abs(E_kin[0]+E_pot[0]), "g", label = "E_total")
    ax.legend()
    plt.ylim([-5.5, 5])
    plt.title('Energy')
    plt.xlabel("t")
    plt.ylabel("E/E0")

    plt.subplot(122)
    plt.plot(Tid, L/np.abs(L[0]))
    plt.title('Angular momentum')
    plt.ylim(-0.1, 1.2)
    plt.xlabel("t")
    plt.ylabel("L/L0")

    plt.show()  # vis figur

def antall_N_pr_runde(vinkel):
    temp = 0
    while temp < N and theta[temp] <= 2*np.pi:
        temp += 1

    return temp


# definere konstanter etc.
k = 1000
m = 1
r0 = np.array([4, 0])
v0 = np.array([0, 10])
t0 = 0
dt = 0.001
N = 10000  # antall iterasjoner

# initsialverdier
r = r0
v = v0
t = t0

X = np.zeros(N)
Y = np.zeros(N)
E_pot = np.zeros(N)  # potensiell energi
E_kin = np.zeros(N)  # kinetisk energi
Tid = np.zeros(N)
L = np.zeros(N)  # angular momentum
theta = np.zeros(N)


# absoluttverdier
r0 = np.linalg.norm(r0)
v0 = np.linalg.norm(v0)

# dimensjonfrie størrelser
R = 1/r0 * r
V = 1/v0 * v
T = (v0 / r0) * t
dT = v0/r0 * dt


for n in range(N):
    # oppdatere verdier
    v_halv = V + (- k/(v0**2 * r0 * m * np.linalg.norm(R)**3))*R * r0/(v0*2)*dT
    R = R + v_halv * r0/v0*dT  # R(t+dt)
    V = v_halv + (- k/(v0**2 * r0 * m * np.linalg.norm(R)**3))* R * r0/(v0*2)*dT  # V(t+dt)
    X[n] = R[0]  # x-koordinat
    Y[n] = R[1]  # y-koordinat
    E_pot[n] = -k*m/(r0 * np.linalg.norm(R))  # potensiell energi
    E_kin[n] = 1/2 * m * v0**2 * np.linalg.norm(V)**2  # kinetisk energi
    Tid[n] = T
    L[n] = m * (R[0] * V[1] - R[1] * V[0])*r0*v0  # angular momentum
    if n == 0:
        theta[n] = 0
    else:  # vinkel
        dTheta = np.abs(np.arccos(X[n-1]/np.sqrt(X[n-1]**2 + Y[n-1]**2)) - np.arccos(X[n]/np.sqrt(X[n]**2 + Y[n]**2)))
        theta[n] = theta[n-1] + dTheta

    T += dT


# one orbit, average energy
N_runde = antall_N_pr_runde(theta)

avg_E_kin = np.sum(E_kin[0:N_runde])/N_runde
avg_E_pot = np.sum(E_pot[0:N_runde])/N_runde

ratio_energy = avg_E_kin/avg_E_pot
print('Forhold mellom kinetisk og potensiell energi: ', ratio_energy)


# sammenligne med eksakt løsning
# r = p/(1 + epsilon * cos(theta))

p = L[0]**2 / (m * k)
epsilon = np.sqrt(1+2*m*(E_kin[0]+E_pot[0])*p/k)

r_eksakt = p/(1+epsilon*np.cos(theta))

X_eksakt = -r_eksakt * np.cos(theta)/r0
Y_eksakt = r_eksakt * np.sin(theta)/r0


# plot figuren
plotter()
