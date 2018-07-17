import numpy as np

class odesolver(object):
    """Different solvers for linear ODEs """
    def __init__(self, start, stop, h, initState, f):
        self.sizeOfState = np.size(initState)
        self.h = h
        self.length = int((stop-start)/self.h) + 1
        self.f = f
        self.initState = initState
        self.t = np.linspace(start, stop, self.length)


    def euler(self):
        """Euler scheme, linear ODEs """
        solution = np.zeros([self.length, self.sizeOfState])
        solution[0, :] = self.initState

        for i in range(self.length - 1):
            solution[i+1, :] = solution[i, :] + self.h * self.f(self.t[i], solution[i, :])

        return solution


    def midpoint(self, h = None):
        """Midpoint scheme, linear ODEs """
        solution = np.zeros([self.length, self.sizeOfState])
        solution[0, :] = self.initState

        for i in range(self.length - 1):
            solution[i+1, :] = solution[i, :] + self.h * self.f(self.t[i] + self.h/2, solution[i, :] + (self.h/2) * self.f(self.t[i], solution[i, :]))

        return solution


    def RK4(self):
        """Runge-Kutta 4 scheme, linear ODEs """
        solution = np.zeros([self.length, self.sizeOfState])
        solution[0, :] = self.initState

        for i in range(self.length - 1):
            k1 = self.f(self.t[i]           , solution[i, :])
            k2 = self.f(self.t[i] + self.h/2, solution[i, :] + (self.h/2) * k1)
            k3 = self.f(self.t[i] + self.h/2, solution[i, :] + (self.h/2) * k2)
            k4 = self.f(self.t[i] + self.h  , solution[i, :] + (self.h) * k3)


            solution[i+1, :] = solution[i, :] + (self.h/6)*(k1 + 2*k2 + 2*k3 + k4)

        return solution
