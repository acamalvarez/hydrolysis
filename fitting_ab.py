import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize
from parameters import Parameters

from functions import coeff_determination, X2, a_function, b_function

reactor = 'small'

a_bData = pd.read_csv('data/a_b.csv')

s0 = a_bData[(a_bData['Reactor'] == reactor)]['s0'] * Parameters.Htot
e0 = a_bData[(a_bData['Reactor'] == reactor)]['e0']

a = a_bData[(a_bData['Reactor'] == reactor)]['a']
b = a_bData[(a_bData['Reactor'] == reactor)]['b']


## fit a

popt, pcov = optimize.curve_fit(a_function, [s0, e0], a, p0=[1, 1, 1],
                bounds=([0,0,0], [1e4, 1e4, 1e4]))

print(popt)

plt.plot([a.min(), a.max()], [a.min(), a.max()], '-')
plt.plot(a, a_function([s0, e0], *popt), 'o', fillstyle='none')
plt.title('a')

r_square = coeff_determination(a, a_function([s0, e0], *popt))
print('R^2 (a) =', r_square)

print('New a:', a_function([s0, e0], *popt))



plt.show()


### fit b

popt, pcov = optimize.curve_fit(b_function, s0, a, p0=[1])


plt.plot([b.min(), b.max()], [b.min(), b.max()], '-')
plt.plot(b, b_function(s0, *popt), 'o', fillstyle='none')
plt.title('b')
plt.show()

r_square = coeff_determination(b, b_function(s0, *popt))
print('R^2 (b) =', r_square)



print('k3 =', popt)

print('New b:', b_function(s0, *popt))
