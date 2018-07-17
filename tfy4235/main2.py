import sys
import os
from odesolver2 import odesolver2
import numpy as np
from scipy.integrate import odeint
from scipy import integrate
import h5py as h5
import time


def options_handeler(args):
    """ Sort the run arguments """

    '''
    if '--help' in args:
        print('\n')
        print('Usage: python3 main2.py [options]\n\n')
        print('    --task:  spesify task (eg. 1.5)\n')
        print('    -h    :  spesify stepsize, default = 0.1\n')
    '''

    if '--task' in args:
        index = args.index('--task')
        task = args[index+1]
    else:
        task = None

    if '-h' in args:
        index = args.index('-h')
        h = float(args[index+1])
    else:
        print("none or invalid input-arguments, spesify with '-h'. default h = 0.1")
        h = 0.1

    return task, h


def Bfield(x, y, z, t):
    """Magnetic field, B """
    B = 1
    beta = 0.2
    return np.array([0, 0, B + beta*y])


def Efield(x, y, z, t):
    """Electric field, E """
    E = 0.1
    return np.array([0, 0, 0])


def f(X, t, dothething = False):
    """ X' = f(X, t) """
    # X = [x1 x2 x3 x4 x5 x6] = [x y z v_x v_y v_z]
    sign = 1
    B = Bfield(X[0], X[1], X[2], t)
    E = Efield(X[0], X[1], X[2], t)
    newstate = np.array([X[3],
                         X[4],
                         X[5],
                         sign*(X[4]*B[2]-X[5]*B[1]),
                         sign*(X[5]*B[0]-X[3]*B[2]),
                         sign*(X[3]*B[1]-X[4]*B[0])])
    newstate += np.array([0, 0, 0, E[0], E[1], E[2]])

    return newstate


def task15(stepsize):
    """ Task 1.5 """

    euler     = True
    midpoint  = True
    RK4       = True

    # initial values
    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    h = stepsize
    X0 = [0, 0, 0, 1, 0, 1]  # initial state [X, V]
    t = np.arange(0, 20, h)

    tic = time.time()
    print('h = {}   N: {}   time: '.format(round(h,4), np.size(t)), end = '\r')
    sol = odesolver2(f, X0, t, h)
    if euler:
        solution_euler = sol.euler()
    if midpoint:
        solution_mid = sol.midpoint()
    if RK4:
        solution_RK4 = sol.RK4()

    print('h = {}   N: {}   time: {}'.format(round(h,4), np.size(t) ,round((time.time()-tic), 5)))

    if euler:
        fil = h5.File('results_euler_h{}.h5'.format(round(h,4)), 'w')
        fil.create_dataset('results', data = solution_euler)
        fil.create_dataset('time', data = t)
        fil.create_dataset('stepsize', data = h)
    if midpoint:
        fil = h5.File('results_midpoint_h{}.h5'.format(round(h,4)), 'w')
        fil.create_dataset('results', data = solution_mid)
        fil.create_dataset('time', data = t)
        fil.create_dataset('stepsize', data = h)
    if RK4:
        fil = h5.File('results_RK4_h{}.h5'.format(round(h,4)), 'w')
        fil.create_dataset('results', data = solution_RK4)
        fil.create_dataset('time', data = t)
        fil.create_dataset('stepsize', data = h)

    fil.close()


def task16(stepsize):
    def exact(t, r_L, omega_c, q, v, pos):
        """ exact solution """
        sol = np.zeros([np.size(t), 3])

        sol[:, 0] = r_L * np.sin(omega_c * t)
        sol[:, 1] = np.sign(q) * r_L * (np.cos(omega_c * t) - 1)
        sol[:, 2] = v[2]*t

        return sol

    def exactReduced(t, r_L, omega_c, q, v, pos):
        """ exact solution reduced units"""
        sol = np.zeros([np.size(t), 3])

        sol[:, 0] = np.sin(t)
        sol[:, 1] = np.sign(q) * (np.cos(t) - 1)
        sol[:, 2] = 1*t

        return sol


    euler     = True
    midpoint  = True
    RK4       = True

    # initial values
    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    B = Bfield(0, 0, 0, 0)
    E = Efield(0, 0, 0, 0)
    q = elementaryCharge
    m = protonMass

    h = stepsize
    X0 = [0, 0, 0, 1, 0, 1]  # initial state [X, V]

    omega_c = abs(q) * np.linalg.norm(B) / m
    r_L = X0[3]/omega_c

    feil_euler = []
    feil_mid = []
    feil_RK4 = []

    h_RK4 = []


    for filename in os.listdir():
        if 'results_' in filename:
            print(filename)
            try:
                f = h5.File(filename, 'r')
                X = f.get('results')
                t = f['time'][()]
                h = f['stepsize'][()]

            except Exception as e:
                raise

            solution_exact = exactReduced(t, r_L, omega_c, q, X0[3:], X0[0:3])
                #calculate the error
            if 'euler' in filename:
                feil_euler.append(np.max(np.linalg.norm(X[:, 0:3], axis = 1) - np.linalg.norm(solution_exact, axis = 1)))


            elif 'midpoint' in filename:
                feil_mid.append(np.max(np.linalg.norm(X[:, 0:3], axis = 1) - np.linalg.norm(solution_exact, axis = 1)))


            elif 'RK4' in filename:
                feil_RK4.append(np.max(np.linalg.norm(X[:, 0:3], axis = 1) - np.linalg.norm(solution_exact, axis = 1)))
                h_RK4.append(h)


    fil = h5.File('error_16.h5', 'w')
    d = fil.create_dataset('error', data = [feil_euler, feil_mid, feil_RK4, h_RK4])
    fil.close()


def task18():
    """ task 1.8 """
    def exactReduced(t, X0):
        """ exact solution reduced units"""
        sol = np.zeros([np.size(t), 3])
        E = 0.1
        B = 1
        v_ExB = E/B

        sol[:, 0] = v_ExB/X0[3] * t + (1 - v_ExB/X0[3])*np.sin(t)
        sol[:, 1] = np.sign(q) * (1 - v_ExB/X0[3])*(np.cos(t) - 1)
        sol[:, 2] = X0[5]/X0[3] *t

        return sol

    euler     = False
    midpoint  = False
    RK4       = True

    # initial values
    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    B = Bfield(0, 0, 0, 0)
    E = Efield(0, 0, 0, 0)
    q = elementaryCharge
    m = protonMass

    h = 0.005
    X0 = [0, 0, 0, 1, 0, 1]  # initial state [X, V]

    omega_c = abs(q) * np.linalg.norm(B) / m
    r_L = X0[3]/omega_c

    t = np.arange(0, 20, h)

    # Exact solution
    solution_exact = exactReduced(t, X0)

    # Numerical solution
    sol = odesolver2(f, X0, t, h)
    solution_RK4 = sol.RK4()

    feil_RK4 = np.max(np.linalg.norm(solution_RK4[:, 0:3], axis = 1) - np.linalg.norm(solution_exact, axis = 1))
    print('Numerical error using RK4: {}'.format(feil_RK4))


    fil = h5.File('task18_results_RK4.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()

    fil = h5.File('task18_results_exa.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()


def task19():
    def kineticEnergy(t, X0, m, q):
        """ Kinetic energy constant E and B """
        K0 = 0.5*(X0[3]**2 + X0[5]**2)
        E = 0.1
        B = 1
        v_ExB = E/B

        K = K0 + np.sign(q)*E*(1-v_ExB/X0[3]) * (np.cos(t)-1)

        return K

    def kineticEnergyNum(t, X, m):
        """ Kinetic energy numerical sol """
        v = X[:, 3::]
        v = np.linalg.norm(v, axis = 1)

        return 0.5 * v**2

    # initial values
    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    B = Bfield(0, 0, 0, 0)
    E = Efield(0, 0, 0, 0)
    q = elementaryCharge
    m = protonMass

    h = 0.1
    X0 = [0, 0, 0, 1, 0, 1]  # initial state [X, V]
    t = np.arange(0, 5, h)

    omega_c = abs(q) * np.linalg.norm(B) / m
    r_L = X0[3]/omega_c

    K = kineticEnergy(t, X0, m, q)

    # Numerical solution
    sol = odesolver2(f, X0, t, h)
    solution_RK4 = sol.RK4()
    K_num = kineticEnergyNum(t, solution_RK4, m)
    solution_euler = sol.euler()
    K_num_euler = kineticEnergyNum(t, solution_euler, m)

    fig = plt.figure()
    plt.plot(t, K, label = 'Exact')
    plt.plot(t, K_num, '.', label = 'RK4')
    plt.plot(t, K_num_euler, 's', alpha = 0.2, label = 'Euler')
    plt.xlabel('t')
    plt.title('Kinetic Energy')
    plt.legend()
    plt.show()


def task112():
    """ Task 1.12 """

    # initial values
    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    h = 0.005
    X0 = [0, 0, 0, 1, 0, 0]  # initial state [X, V]
    t = np.arange(0, 40, h)

    sol = odesolver2(f, X0, t, h)
    solution_RK4 = sol.RK4()

    '''
    Man må manulelt gå inn i f() å endre fortegn
    '''

    fil = h5.File('task112_results_RK4_pos.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()


def task114():
    """Task 1.14 """

    def Bfield(x, y, z, t):
        """Magnetic field, B """
        B = 1
        r = np.sqrt(x**2 + y**2)
        if abs(x) < 0.0001 and abs(y) < 0.0001:
            phi = 0
        elif x >= 0.0001:
            phi = np.arcsin(y/r)
        elif x > 0.0001:
            phi = np.arctan(y/x)
        elif x < -0.0001:
            phi = -np.arcsin(y/r) + np.pi


        return B * np.array([-np.sin(phi), np.cos(phi), 0])

    def Bfield2(x, y, z, t):
        """Magnetic field, B """
        B = 5
        r = np.sqrt(x**2 + y**2)
        theta = np.arccos(z/r)
        if abs(x) < 0.001 and y >= 0:
            phi = np.pi/2
        elif abs(x) < 0.001 and y < 0:
            phi = -np.pi/2
        else:
            phi = np.arctan(y/x)
        if abs(phi) > 10:
            print(phi)

        return B * np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)])


    def Efield(x, y, z, t):
        """Electric field, E """
        E = 1
        return np.array([0, 0, 0])


    def f(X, t):
        """ X' = f(X, t) """
        # X = [x1 x2 x3 x4 x5 x6] = [x y z v_x v_y v_z]
        sign = 1
        B = Bfield(X[0], X[1], X[2], t)
        E = Efield(X[0], X[1], X[2], t)
        newstate = np.array([X[3],
                             X[4],
                             X[5],
                             sign*(X[4]*B[2]-X[5]*B[1]),
                             sign*(X[5]*B[0]-X[3]*B[2]),
                             sign*(X[3]*B[1]-X[4]*B[0])])
        newstate += np.array([0, 0, 0, E[0], E[1], E[2]])

        return newstate


    h = 0.001
    X0 = [6, 0, 0, 0, np.sqrt(2)/2, np.sqrt(2)/2]  # initial state [X, V]
    t = np.arange(0, 90, h)

    # sol = odesolver2(f, X0, t, h)
    # solution_RK4 = sol.RK4()
    solution_RK4 = odeint(f, X0, t)
    '''
    Man må manulelt gå inn i f() å endre fortegn
    '''

    fil = h5.File('task114_results_RK4_pos.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()


def hem2():
    """helmholtz"""
    def Bfield(x, y, z, t):
        """Magnetic field, B """

        r = (x*x + y*y)**(0.5)

        def br(r, z):
            y1 = lambda theta: (z-1)*np.cos(theta)/(( (r-R*np.cos(theta))**2 + R**2 * np.sin(theta)**2 + (z-1)**2 )**(3/2))

            y2 = lambda theta: (z+1)*np.cos(theta)/(( (r-R*np.cos(theta))**2 + R**2 * np.sin(theta)**2 + (z+1)**2 )**(3/2))

            I1 = integrate.quad(y1, 0, 2*np.pi)
            I2 = integrate.quad(y2, 0, 2*np.pi)

            return(I1 + I2)

        def bz(r, z):
            y1 = lambda theta: (1-(r/R)*np.cos(theta))/(( (r-R*np.cos(theta))**2 + R**2 * np.sin(theta)**2 + (z-1)**2 )**(3/2))

            y2 = lambda theta: (1-(r/R)*np.cos(theta))/(( (r-R*np.cos(theta))**2 + R**2 * np.sin(theta)**2 + (z+1)**2 )**(3/2))

            I1 = integrate.quad(y1, 0, 2*np.pi)
            I2 = integrate.quad(y2, 0, 2*np.pi)

            return(I1 + I2)

        B_r = ((1+R**2)**(3/2))/(4*np.pi * R) * br(r, z)[0]
        B_z = ((1+R**2)**(3/2))/(4*np.pi) * bz(r, z)[0]

        return (1/B0) * np.array([B_r * x/r, B_r*y/r, B_z])


    def f(X, t, dothething = False):
        """ X' = f(X, t) """
        # X = [x1 x2 x3 x4 x5 x6] = [x y z v_x v_y v_z]
        sign = 1
        B = Bfield(X[0], X[1], X[2], t)
        if dothething:
            normB.append(np.linalg.norm(B))
            kinetic_energy.append(X[3]**2 + X[4]**2 + X[5]**2)
        #E = Efield(X[0], X[1], X[2], t)
        newstate = np.array([X[3],
                             X[4],
                             X[5],
                             (X[4]*B[2]-X[5]*B[1]),
                             (X[5]*B[0]-X[3]*B[2]),
                             (X[3]*B[1]-X[4]*B[0])])
        #newstate += np.array([0, 0, 0, E[0], E[1], E[2]])
        print('time: {}/{}'.format(round(t, 1), totalTime), end='\r')

        if dothething:
            vnorm = X[3::] - np.inner(X[3::], B)*B/(normB[-1]**2)
            mu.append(-(np.linalg.norm(vnorm)**2)/normB[-1])

        return newstate

    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    h = 0.001
    B0 = 0.001
    d = 1
    R = 0.2

    totalTime = 50

    omega_c = elementaryCharge * B0 / electronMass


    X0 = [1, 0, 0, np.sqrt(2)/2, 0, np.sqrt(2)/2]
    normB = []
    mu = []
    kinetic_energy = []
    # X0 = [6, 0, 0, 0, np.sqrt(2)/2, np.sqrt(2)/2]  # initial state [X, V]
    #X0 = [6, 0, 0, 1, 0, 0]  # initial state [X, V]
    t = np.arange(0, totalTime, h)

    sol = odesolver2(f, X0, t, h)
    solution_RK4 = sol.RK4()
    #solution_RK4 = odeint(f, X0, t)
    print('time: {}'.format(time.time()-tic))

    fil = h5.File('task21Hem3_results_RK4_v2.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()

    fil = h5.File('task21Hem3_magnetic_moment.h5', 'w')
    fil.create_dataset('mu', data = mu)
    fil.create_dataset('time', data = t)
    fil.create_dataset('kinetic_energy', data = kinetic_energy)
    fil.close()


def earth():
    """earth's radiationbelt"""
    def Bfield(x, y, z, t):
        """B-field"""

        r = (x*x + y*y + z*z)**(0.5)
        Bx = -3*x*z/(r**5)
        By = -3*y*z/(r**5)
        Bz = -(1/(r**5)) * (2*z*z - x*x -y*y)

        return np.array([Bx, By, Bz])

    def f(X, t, dothething=False):
        """ X' = f(X, t) """
        # X = [x1 x2 x3 x4 x5 x6] = [x y z v_x v_y v_z]
        # sign = 1
        B = Bfield(X[0], X[1], X[2], t)
        if dothething:
            normB.append(np.linalg.norm(B))
            kinetic_energy.append(X[3]**2 + X[4]**2 + X[5]**2)
        newstate = np.array([X[3],
                             X[4],
                             X[5],
                             alpha*(X[4]*B[2]-X[5]*B[1]),
                             alpha*(X[5]*B[0]-X[3]*B[2]),
                             alpha*(X[3]*B[1]-X[4]*B[0])])
        # newstate += np.array([0, 0, 0, E[0], E[1], E[2]])
        print('time: {}/{}'.format(round(t, 1), totalTime), end='\r')

        if dothething:
            vnorm = X[3::] - np.inner(X[3::], B)*B/(normB[-1]**2)
            mu.append(-(np.linalg.norm(vnorm)**2)/normB[-1])

        return newstate


    radiusEarh = 6378000
    B0 = 30.7e-6
    c0 = 299792458

    protonMass          = np.float(1.673e-27)
    electronMass        = np.float(9.109e-31)
    elementaryCharge    = np.float(1.602e-19)

    omega = c0/radiusEarh
    omega2 = B0 * elementaryCharge/protonMass

    alpha = omega2/omega

    normB = []
    mu = []
    kinetic_energy = []

    #h = 0.01
    h = 0.1
    totalTime = 10000
    X0 = [0, 4, 0, 0.08*np.sqrt(2)/2, 0, 0.08*np.sqrt(2)/2]
    # X0 = [6, 0, 0, 0, np.sqrt(2)/2, np.sqrt(2)/2]  # initial state [X, V]
    #X0 = [6, 0, 0, 1, 0, 0]  # initial state [X, V]
    t = np.arange(0, totalTime, h)

    tic = time.time()
    sol = odesolver2(f, X0, t, h)
    solution_RK4 = sol.RK4()
    #solution_RK4 = odeint(f, X0, t)
    print('time: {}'.format(time.time()-tic))

    fil = h5.File('task3_results_RK4_v2.h5', 'w')
    fil.create_dataset('results', data = solution_RK4)
    fil.create_dataset('time', data = t)
    fil.create_dataset('stepsize', data = h)
    fil.close()

    fil = h5.File('task3_magnetic_moment.h5', 'w')
    fil.create_dataset('mu', data = mu)
    fil.create_dataset('time', data = t)
    fil.create_dataset('kinetic_energy', data = kinetic_energy)
    fil.close()



def main():

    task, h = options_handeler(sys.argv)

    if task is None:
        print("No task")
        # options_handeler('--help')

    elif task == '1.5':
        print('\nrunning task 1.5')
        task15(h)

    elif task == '1.6':
        print('\nrunning task 1.6')
        task16(h)

    elif task == '1.8':
        print('\nrunning task 1.8')
        task18()

    elif task == '1.9':
        print('\nrunning task 1.9')
        task19()

    elif task == '1.12':
        print('\nrunning task 1.12')
        task112()

    elif task == '1.14':
        print('\nrunning task 1.14')
        task114()

    elif task == 'hem2':
        print('\nrunning task hem')
        hem2()

    elif task == '3':
        print('\nrunning task 3')
        earth()


if __name__ == '__main__':
    main()
