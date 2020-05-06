import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

from Hydrolysis_functions import Qi, Fit_Hydrolysis, Hydrolysis_file, Files

x = Hydrolysis_file(Files.Large.s10e3).time_sec
y = Hydrolysis_file(Files.Large.s10e3).Volume

print(x)

f = interp1d(x, y)

x_new = np.linspace(0, 11760, num=int(11760/5+1))
y_new = f(x_new)

print(x_new)


plt.plot(x, y, 'o', fillstyle='none')
# plt.plot(x_new, y_new+2, 'o', fillstyle='none')

plt.show()