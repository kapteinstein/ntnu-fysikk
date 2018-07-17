import matplotlib.pyplot as plt
import numpy as np

m = np.array([0,0,1])
boks = 5
N = 60

x, z = np.meshgrid(np.linspace(-boks, boks, N), np.linspace(-boks, boks, N))

B_x = (3 * x * (m[0] * x + m[2] * z))/((x**2 + z**2)**(5/2)) - m[0]/((x**2 + z**2)**(3/2))
B_z = (3 * z * (m[0] * x + m[2] * z))/((x**2 + z**2)**(5/2)) - m[2]/((x**2 + z**2)**(3/2))

length  = np.sqrt(x**2 + z**2)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-boks, boks), ylim=(-boks, boks))

ax.scatter(0, 0, s=1000, color='k')
ax.quiver(x, z, B_x, B_z)

plt.show()