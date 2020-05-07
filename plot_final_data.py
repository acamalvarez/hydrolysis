from functions import plot_final
from parameters import Parameters
from files import Files
import numpy as np



files = [Files.Small.s5e1, Files.Small.s10e1, Files.Small.s15e1, Files.Small.s30e1]
s0s = [5 * Parameters.Htot,  10 * Parameters.Htot,  15 * Parameters.Htot,  30 * Parameters.Htot]
legends = ['5 g/L', '10 g/L', '15 g/L', '30 g/L']


plot_final(files, s0s, 1, legends)

