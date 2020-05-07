import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0, 100, 0.1)

a_values = [2.234, 2.4]
b_values = [0.789, 0.9]




for i in np.arange(len(a_values)):
    P = (1 / b_values[i]) * np.log(a_values[i] * b_values[i] * t + 1)

    plt.plot(t, P)

plt.show()