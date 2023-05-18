import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
from os.path import *


file_loc = sys.argv[1]

times = []
g2 = []
with open(file_loc) as file:
    for line in file:
        if "Computing second order cumulant lineshape function" in line:
            line = next(file)
            line = next(file)
            while True > 0:
                print("Line: ", line[0:-1])
                line = next(file)
                sp = line.split()
                if len(sp) == 0:
                    break
                times.append(float(sp[1]))
                g2.append(float(sp[2]))
            break
times = np.array(times)
g2 = np.array(g2)

fig, ax = plt.subplots()
ax.plot(times, g2)

ax.set_xlabel('Time (fs)')
ax.set_ylabel('Re[g_2]')
fig.tight_layout()
plt.show()
                