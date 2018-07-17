import numpy as np

class odesolver2(object):
    """odesolver MK-2"""
    def __init__(self, f, state, t, h):
        self.f = f
        self.state = state
        self.t = t
        self.h = h

    def euler(self):
        """Euler metode"""
        solution = np.zeros([np.size(self.t), np.size(self.state)])
        solution[0, :] = self.state

        for i in range(np.size(self.t)-1):
            solution[i + 1, :] = solution[i, :] + self.h * self.f(solution[i, :], self.t[i])

        return solution


    def midpoint(self):
        """Midpoint scheme"""
        solution = np.zeros([np.size(self.t), np.size(self.state)])
        solution[0, :] = self.state

        for i in range(np.size(self.t)-1):
            solution[i+1, :] = solution[i, :] + self.h * self.f(solution[i, :] + (self.h/2) * self.f(solution[i, :], self.t[i]), self.t[i] + self.h/2)

        return solution

    def RK4(self):
        """Runge-Kutta scheme"""
        solution = np.zeros([np.size(self.t), np.size(self.state)])
        solution[0, :] = self.state

        for i in range(np.size(self.t)-1):
            k1 = self.f(solution[i, :]                  , self.t[i])
            k2 = self.f(solution[i, :] + (self.h/2) * k1, self.t[i] + self.h/2)
            k3 = self.f(solution[i, :] + (self.h/2) * k2, self.t[i] + self.h/2)
            k4 = self.f(solution[i, :] + (self.h) * k3  , self.t[i] + self.h)

            self.f(solution[i, :], self.t[i], dothething = True)


            solution[i+1, :] = solution[i, :] + (self.h/6)*(k1 + 2*k2 + 2*k3 + k4)

        return solution
