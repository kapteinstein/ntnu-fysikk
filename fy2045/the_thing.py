# Numerisk prosjekt, Kvantefysikk 1 (FY2045), NTNU 2015

import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation



'''
# wave function = wave packet:
# psi = Ce^(-(x-x0)^2/(2*sigma_x**2))*e^i(kx-wt)
#     = Ce^(-(x-x0)^2/(2*sigma_x**2))*(cos(kx-wt) + i*sin(kx-wt))    # hva er C ??
#
# psi = psi_r + i * psi_i
# psi_r = Ce^(-(x-x0)^2/(2*sigma_x**2)) * ( cos(kx-wt) )
# psi_i = Ce^(-(x-x0)^2/(2*sigma_x**2)) * ( sin(kx-wt) )

'''

def Psi(x, x0, w, C, sigma_x, dt):
    psi_r = C * np.exp(-(x - x0)**2 / (2*sigma_x**2)) * (np.cos(k*x))
    psi_i = C * np.exp(-(x - x0)**2 / (2*sigma_x**2)) * (np.sin(k*x - w*dt/2))
    psi_r[0] = 0
    psi_r[-1] = 0
    psi_i[0] = 0
    psi_i[-1] = 0
    return psi_r + 1j*psi_i


# bølgefunksjonen på et senere tidspunkt, ligning 7a,b.
def Psi_dt_v2(a, dt, number_of_times):
    Nx = np.size(a)
    for i in range(number_of_times):
        a.imag[1:Nx-1] = a.imag[1:Nx-1] + dt * (a.real[1:Nx-1] * V[1:Nx-1]/hbar
                            + hbar*(a.real[2:Nx] - 2*a.real[1:Nx-1] + a.real[0:Nx-2]) / (2*m*(dx**2)))
        a.real[1:Nx-1] = a.real[1:Nx-1] - dt * (a.imag[1:Nx-1] * V[1:Nx-1]/hbar
                            + hbar*(a.imag[2:Nx] - 2*a.imag[1:Nx-1] + a.imag[0:Nx-2]) / (2*m*(dx**2)))
    return a

# definer potensialet
def potential(a, level):
    b = np.zeros(np.size(a))
    for n in range(np.size(a) - 1):
        if (L/2 - l/2) < a[n] < (L/2 + l/2):
            b[n] = level * E
    return b

def newton_integral(psi):
    prob = (psi.real**2 + psi.imag**2)
    integral1 = 0
    integral2 = 0
    for i in range(int(N/2) - 1):
        integral1 += prob[i + 1] + prob[i]
        integral2 += prob[int(N/2) + i + 1] + prob[int(N/2) + i]
    return integral1 * 0.5*L/N, integral2 * 0.5*L/N


def plotter(X, a):  # definer plottefunksjonen.
    plt.clf()  # sletter forrige graf
    plt.plot(X, a.real, "b", label = "$\Psi_R$")
    plt.plot(X, a.imag, "g", label = "$\Psi_I$")  # plotter psi
    plt.plot(X, (a.real**2 + a.imag**2), "r",  label = "$|\Psi|^2$")
    plt.plot(X, V/E, "k", linewidth = 2, label = "$V$")  # plotter potensialet
    plt.legend()
    plt.xlim((0, 20))  # setter aksene
    plt.ylim((-2, 2))
    plt.pause(0.001)  # trenger denne for 'animasjon' kan kanskje tas bort hvis kun ett plot.



# definer variabler, funksjoner, etc.
pi = np.pi  # gode gamle pi = 3.1415...
hbar = 1  # h-bar
m = 1  # mass, m
k = 20  # wavenumber
L = 20  #
sigma_x = 1
C = np.sqrt(np.sqrt(1/(pi*sigma_x**2)))  # normalization constant
w = (k**2/2)  # omega
l = L/50  # lengde på barriere
E = hbar*w  # energy


# definer området denne shiten skal virke over, eg. fra x0 -> x1
N = 1000  # antall punkter på x-aksen
x0 = 5  # startpunkt for psi
x1 = L  # slutt
t = 0  # starttid, overflødig
dx = L/(N-1)
dt = 0.01 * 2*m*(dx**2)  # for stabilitet
T = int(L*m/(2*k*dt))


oppgave_nr = 50


###########   OPPGAVE 1   #############
if oppgave_nr == 1:
    # initsialisere, dvs psi(x0, t = 0)
    X = np.linspace(0, x1, num=N)
    V = potential(X, 0)
    psi = Psi(X, x0, w, C, sigma_x, dt)  # går sikkert ann å gjøre denne penere eller på en bedre måte.

    plotter(X, psi)
    plt.savefig('oppgave_nr1.eps')
    plt.show()  # vise figuren


###########   OPPGAVE 2   #############
if oppgave_nr == 2:
    # initsialisere, dvs psi(x0, t = 0)
    X = np.linspace(0, x1, num=N)
    V = potential(X, 0)
    psi = Psi(X, x0, w, C, sigma_x, dt)  # går sikkert ann å gjøre denne penere eller på en bedre måte.

    plt.figure()
    plt.plot(X, psi.real**2+psi.imag**2, 'r', label = "$|\Psi_0|^2$")
    psi_ny = Psi_dt_v2(psi, dt, T)
    plt.plot(X, psi_ny.real**2 + psi_ny.imag**2, 'b', label = "$|\Psi|^2$")
    plt.legend()
    #plt.savefig('oppgave2_1.pdf')
    plt.show()

    '''
    # kjør animasjon
    plt.ion() # må ha med denne for å kunne modifisere plottet etter at det er tegnet opp.
    for n in range(T):
        psi = Psi_dt_v2(psi, dt, T)
        plotter(X, psi)
        plt.show()  # vise figuren
    '''

###########   OPPGAVE 3   #############
elif oppgave_nr == 3:
    # initsialisere, dvs psi(x0, t = 0)
    X = np.linspace(0, x1, num=N)
    V = potential(X, 0.5)
    psi = Psi(X, x0, w, C, sigma_x, dt)  # går sikkert ann å gjøre denne penere eller på en bedre måte.

    psi = Psi_dt_v2(psi, dt, T)

    Reflection, Transmission = newton_integral(psi)

    print('R: ', Reflection)
    print('T: ', Transmission)
    print('R + T: ', Reflection + Transmission)

    plotter(X, psi)
    plt.savefig('oppgave3.pdf')
    plt.show()  # vise figuren

###########   OPPGAVE 4   #############
elif oppgave_nr == 4:

    X = np.linspace(0, x1, num=N)

    Reflection = np.zeros(50)
    Transmission = np.zeros(50)

    for i in range(50):
        V = potential(X, 3*i/100)
        psi = Psi(X, x0, w, C, sigma_x, dt)
        psi = Psi_dt_v2(psi, dt, T)
        Reflection[i], Transmission[i] = newton_integral(psi)
        psi = None
        print('%i of 50' % (i+1))

    plt.plot(np.linspace(0, 3/2, 50), Reflection, "r", label = "R")
    plt.plot(np.linspace(0, 3/2, 50), Transmission, "b", label = "T")
    plt.xlabel("$V_0/E$")
    plt.legend()
    plt.savefig('oppgave4.pdf')
    plt.show()  # vise figuren


###########   OPPGAVE 5   #############
elif oppgave_nr == 5:
    #N = 2000
    X = np.linspace(0, x1, num=N)

    Reflection = np.zeros(50)
    Transmission = np.zeros(50)
    #plt.ion()

    for i in range(50):
        l = i*L/1000
        V = potential(X, 9/10)
        psi = Psi(X, x0, w, C, sigma_x, dt)
        psi = Psi_dt_v2(psi, dt, T)
        #plotter(X, psi)
        #plt.show()  # vise figuren
        Reflection[i], Transmission[i] = newton_integral(psi)
        print('%i of 50' % (i+1))


    plt.plot(np.linspace(0, L/20, 50), Reflection, np.linspace(0, L/20, 50), Transmission)
    plt.savefig('oppgave5.pdf')
    plt.show()  # vise figuren


else:
    print('\n https://www.youtube.com/watch?v=Znby3t3AS5s')