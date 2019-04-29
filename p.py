

import numpy as np
import matplotlib.pyplot as plt

fsolution=np.loadtxt('mat.txt',delimiter=",")
xmin = 0.0
xmax = 1.0

ymin = 0.0
ymax = 1.0
plt.imshow(np.transpose(fsolution),
                   origin="lower", interpolation="nearest",
                   extent=[xmin, xmax, ymin, ymax])

plt.xlabel("x")
plt.ylabel("y")

plt.colorbar()

plt.tight_layout()
plt.savefig("t.png")
