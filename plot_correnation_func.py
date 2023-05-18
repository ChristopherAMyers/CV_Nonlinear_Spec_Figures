import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
from os.path import *


file_loc = sys.argv[1]
dt = 4.0
g2 = np.loadtxt(file_loc)
times = np.arange(len(g2))*dt

fig, ax = plt.subplots()
ax.plot(times, g2)

ax.set_xlabel('Time (fs)')
ax.set_ylabel('Re[g_2]')
fig.tight_layout()
plt.show()
                