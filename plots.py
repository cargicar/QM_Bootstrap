import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
size, grid_size, g =np.loadtxt('fort.1')
x2_energy=np.transpose(np.loadtxt('fort.2'))


fig, ax = plt.subplots()

ax.scatter(x2_energy[0],x2_energy[1])
ax.set_title(f"Anharmonic  k={size}, g={g}, grid size = {grid_size}.png")
#ax.set(xlim=(1.36, 1.41), ylim=(0.29, 0.35))
ax.set_xlabel('Energy')
ax.set_ylabel('$< x^2 > $')
plt.savefig(f"plots/fortran Data, k={size}, g={g}, grid size = {grid_size}.png")

