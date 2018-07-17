import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# vector field

N = 8
m = np.array([0, 0, 1])
mu = 100

x, y, z = np.meshgrid(np.linspace(-3, 3, N), np.linspace(-3, 3, N), np.linspace(-3, 3, N))  # startpunkt

r = np.array([x, y, z])
B = r
for i in range(3):
    B[i] = (3 * r[i] * (m[0] * x + m[1] * y + m[2] * z))/((x**2 + y**2 + z**2)**(5/2))\
           - m[i]/((x**2 + y**2 + z**2)**(3/2))



# plotte klode + magnetfelt
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')

theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)
R = 1

X_kule = R * np.outer(np.cos(phi), np.sin(theta))
Y_kule = R * np.outer(np.sin(phi), np.sin(theta))
Z_kule = R * np.outer(np.ones(np.size(theta)), np.cos(theta))

ax1.streamplot(x, y, z, B[0], B[1], B[2])
ax1.plot_surface(X_kule, Y_kule, Z_kule, rstride=2, cstride=2, color='w')


# 2D plot
# se andre fil

plt.show()