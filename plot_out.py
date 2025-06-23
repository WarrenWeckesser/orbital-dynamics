import numpy as np
import matplotlib.pyplot as plt



w = np.loadtxt('out')

plt.figure(1)
plt.clf()
for i in range(len(w)-34,  len(w)):
    plt.plot(w[i, [1, 3, 5, 1]], w[i, [2, 4, 6, 2]], 'k', alpha=0.1)
    plt.plot(w[i, 1], w[i, 2], 'ro')
    plt.plot(w[i, 3], w[i, 4], 'go')
    plt.plot(w[i, 5], w[i, 6], 'bo')

plt.plot(0, 0, 'ko', ms=8)
plt.grid()
plt.axis('equal')

plt.show()
