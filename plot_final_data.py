from functions import plot_final
from parameters import Parameters
from files import Files




### s0 = 5, e0 = 1
plot_final(Files.Small.s5e1, 5 * Parameters.Htot, 1)

### s0 = 10, e0 = 1
plot_final(Files.Small.s10e1, 10 * Parameters.Htot, 1)

### s0 = 15, e0 = 1
plot_final(Files.Small.s15e1, 15 * Parameters.Htot, 1)

### s0 = 30, e0 = 1
plot_final(Files.Small.s30e1, 30 * Parameters.Htot, 1)

### s0 = 40, e0 = 1
plot_final(Files.Small.s40e1, 40 * Parameters.Htot, 1)