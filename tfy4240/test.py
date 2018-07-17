# Numerisk prosjekt, TFY4240, Erik Liodden, Oppgave 1

import numpy as np
import matplotlib
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import time


def potential(x, y, i, a_n):  # potential
    v = np.zeros([np.size(x),np.size(y)])
    for m in range(np.size(y)):
        for n in range(np.size(x)):
            for j in range(i):
                v[n, m] += a_n[j] * np.sin((j + 1) * np.pi * x[m]) \
                           * np.sinh((j + 1) * np.pi * y[n])
    return v


def a_v4(antall, nr):  # fourier constant
    temp = np.linspace(0, 1, N_simpson)
    h = 1/N_simpson
    odde_ledd = 0
    like_ledd = 0
    a = np.zeros(antall)

    for n in range(antall):
        n += 1
        for j in range(int(N_simpson/2 - 1)):
            j += 1
            like_ledd += v0(temp[2*j], nr) * np.sin(n*np.pi*temp[2*j])

        for j in range(int(N_simpson/2)):
            j += 1
            odde_ledd += v0(temp[2*j-1], nr) * np.sin(n*np.pi*temp[2*j-1])

        a[n-1] = (2/(np.sinh(n*np.pi)))* h/3 * (v0(temp[0], nr)
                    * np.sin(n*np.pi*temp[0]) + 2*like_ledd
                        + 4*odde_ledd + v0(temp[-1], nr)
                                                  * np.sin(n*np.pi*temp[-1]))
        like_ledd = 0
        odde_ledd = 0

    return a


def v0(x_kordinater, hvilket_potensial):  # initial potential
    if hvilket_potensial == 1:
        return np.sin(3*np.pi*x_kordinater)
    elif hvilket_potensial == 2:
        return (1-(x_kordinater - 1/2)**4)
    elif hvilket_potensial == 3:
        if x_kordinater -1/2 < 0 or 3/4 - x_kordinater < 0:
            return 0
        else:
            return 1
    else:
        return 1


def electric_field(x, y, a_n):  # electric field
    u = a_n[0] * np.pi * np.cos(np.pi*x) * np.sinh(np.pi*y)
    v = a_n[0] * np.pi * np.sin(np.pi*x) * np.cosh(np.pi*y)
    for i in range(antall_ledd -1):
        i = i+1
        u = u + a_n[i] * (i+1)*np.pi * np.cos((i+1)*np.pi*x) \
                * np.sinh((i+1)*np.pi*y)
        v = v + a_n[i] * (i+1)*np.pi * np.sin((i+1)*np.pi*x) \
                * np.cosh((i+1)*np.pi*y)
    return -u, -v


###########################
#                         #
#  Definer variabler etc  #

###########################

antall_ledd = 50  # hvor mange ledd man skal ta med av V og simpson
N_simpson = 1000  # antall ledd i simpsons metode
N = 100  # hvor mange punkter man skal ha med

'''
potensial nr.1:
V = sin(npix)

potensial nr 2:
v = 1-(x-1/2)**4

potensial nr 3:
v = step-potential between 1/2 and 3/4.
'''

# specify which potential to use
potential_nr = 2
save_figures = False

# define the box
X = np.linspace(0, 1, N)
Y = np.linspace(0, 1, N)

# calculate the potential inside the box
a_n = a_v4(antall_ledd, potential_nr)
pot = potential(X, Y, antall_ledd, a_n)

# calculate the electric field
X1, Y1 = np.meshgrid(np.linspace(0, 0.98, 20), np.linspace(0, 0.98, 20))
X11 = np.linspace(0, 0.98, 20)
Y11 = np.linspace(0, 0.98, 20)
pot2 = potential(X11, Y11, antall_ledd, a_n)
E_x, E_y = electric_field(X1, Y1, a_n)

# normalize the field vectors
E_norm = np.sqrt(E_x**2 + E_y**2)
E_norm[E_norm == 0] = 1

E_x_norm = E_x/E_norm
E_y_norm = E_y/E_norm


# plot potential
fig = plt.figure()
levels = np.linspace(np.min(pot), np.max(pot), 1000)
plt.contourf(X, Y, pot, levels=levels)
plt.colorbar()
plt.title('Potential, $V/V_c$')
plt.xlabel(r"$\xi$")
plt.ylabel("$\eta$")
if save_figures:
    plt.savefig('potential_%i.pdf' % potential_nr)  # save figure

# plot electric field
fig = plt.figure()
levels2 = np.linspace(-np.max(pot), np.max(pot), 1000)
plt.quiver(X1, Y1, E_x_norm, E_y_norm, pot2, cmap=cm.seismic, headlength=7)
plt.colorbar()
plt.title(r'Electric field, $\vec{E}$')
plt.xlabel(r"$\xi$")
plt.ylabel("$\eta$")
if save_figures:
    plt.savefig('electricField_%i.pdf' % potential_nr)  # save figure


# plot potensialer på kanten

plt.figure()
plt.subplots_adjust(hspace=0.4, wspace=0.5)  # spacing between the subplots

# y = 0
plt.subplot(221)
plt.plot(X, pot[0,:])
plt.ylim([-3e-16, 3e-16])
plt.title(r'1: $V(\xi,\,\eta = 0)$')
plt.xlabel(r"$\xi$")
plt.ylabel(r"$V(\xi)$")

# y = L
plt.subplot(222)
plt.plot(X, pot[-1,:])
if potential_nr == 3:
    plt.plot(np.linspace(1/2, 3/4, 100), np.ones(100))
else:
    plt.plot(X, v0(X,potential_nr))
plt.title(r'2: $V(\xi,\,\eta = 1)$')
plt.xlabel(r"$\xi$")
plt.ylabel(r"$V(\xi)$")

# x = 0
plt.subplot(223)
plt.plot(X, pot[:,0])
plt.title(r'3: $V(\xi = 0,\,\eta)$')
plt.xlabel("$\eta$")
plt.ylabel("$V(\eta)$")

# x = L
plt.subplot(224)
plt.plot(X, pot[:,-1])
plt.title(r'4: $V(\xi = 1,\,\eta)$')
plt.xlabel("$\eta$")
plt.ylabel("$V(\eta)$")

if save_figures:
    plt.savefig('boundary_%i.pdf' % potential_nr)  # save figure


# calculate error
# I know it takes forever, but it gets the job done.

diff = np.zeros(antall_ledd-1)
for i in range(1, antall_ledd-1):
    pot1err = potential(X, Y, i, a_n)
    pot2err = potential(X, Y, i+1, a_n)
    diff[i] = np.sum((pot1err - pot2err)**2)
    print('%i av %i' % (i+2, antall_ledd))
    pot1err = None
    pot2err = None
    b_n = None

plt.figure()
plt.plot(range(antall_ledd-1), diff/np.max(diff), '.')
plt.yscale('log')
plt.title('Feil bellom ett ledd og det neste')
plt.ylabel('Relativ feil')
plt.xlabel('Antall Fourier ledd i potensialet ')
if save_figures:
    plt.savefig('feil_%i.pdf' % potential_nr)  # save figure

plt.show()
