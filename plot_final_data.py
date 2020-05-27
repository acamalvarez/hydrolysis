from functions import plot_final
from Qi_functions import fit_several
from parameters import Parameters
from files import Files
import numpy as np


### Uncomment to see the final graph with one parameter for all the hydrolysis curves

files = [Files.Small.s5e1, Files.Small.s10e1, Files.Small.s15e1]
s0s = [5 * Parameters.Htot,  10 * Parameters.Htot,  15 * Parameters.Htot]

legends = ['45.5 mM', '91.0 mM', '136.5 mM']

plot_final(files, s0s, 1, legends)

# files = [Files.Small.s5e1, Files.Small.s10e1, Files.Small.s15e1, Files.Small.s30e1]

# s0s = np.array([5, 10, 15, 30]) * Parameters.Htot
# e0s = np.array([1, 1, 1, 1, 1])

# fit_several(files, function='p_Valencia', s0s=s0s, e0s=e0s, reactor='small')

