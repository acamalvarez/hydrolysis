import numpy as np
import matplotlib.pyplot as plt

from functions import plot_several
from Qi_functions import Qi, fit_Qi
from our_functions import Function, fit_function
from process_file import Hydrolysis_file
from files import Files
from parameters import Parameters

tC2NI = np.array([0, 1500, 3000, 4500, 6000, 7500, 9000, 10500, 12000]) / 60
hC2NI = np.array([0, 5.723011831, 7.09950903, 8.130385719,
                  8.775664921, 9.299014378, 9.695530863, 9.957369036, 10.21381354])
# fit_Qi(t=tC2NI, h=hC2NI, function='no', e0=2.96, s0=30)
# fit_Qi(t=tC2NI, h=hC2NI, function='sp', e0=2.96, s0=30)


# plot_several([Files.Small.s5e1, Files.Small.s10e1, Files.Small.s15e1,
#             Files.Small.s20e1, Files.Small.s40e1], title='',
#             legend=['5 g/L', '10 g/L', '15 g/L', '20 g/L', '40 g/L'])



# fit_Qi(t=Hydrolysis_file(Files.Small.s5e1).time_min,
#                P=Hydrolysis_file(Files.Small.s5e1).alphaNH(),
#                function='ab_Valencia', s0=5*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s10e1).time_min,
#                P=Hydrolysis_file(Files.Small.s10e1).alphaNH(),
#                function='ab_Valencia', s0=10*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s15e1).time_min,
#                P=Hydrolysis_file(Files.Small.s15e1).alphaNH(),
#                function='ab_Valencia', s0=15*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s20e1).time_min,
#                P=Hydrolysis_file(Files.Small.s20e1).alphaNH(),
#                function='ab_Valencia', s0=20*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s25e1).time_min,
#                P=Hydrolysis_file(Files.Small.s25e1).alphaNH(),
#                function='ab_Valencia', s0=25*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s30e1).time_min,
#                 P=Hydrolysis_file(Files.Small.s30e1).alphaNH(),
#                 function='ab_Valencia', s0=30*Parameters.Htot, e0=1)

# fit_Qi(t=Hydrolysis_file(Files.Small.s40e1).time_min,
#                 P=Hydrolysis_file(Files.Small.s40e1).alphaNH(),
#                 function='ab_Valencia', s0=40*Parameters.Htot, e0=1)





fit_Qi(t=Hydrolysis_file(Files.Small.s5e1).time_min,
               P=Hydrolysis_file(Files.Small.s5e1).alphaNH(),
               function='p_Valencia', s0=5*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s10e1).time_min,
               P=Hydrolysis_file(Files.Small.s10e1).alphaNH(),
               function='p_Valencia', s0=10*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s15e1).time_min,
               P=Hydrolysis_file(Files.Small.s15e1).alphaNH(),
               function='p_Valencia', s0=15*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s20e1).time_min,
               P=Hydrolysis_file(Files.Small.s20e1).alphaNH(),
               function='p_Valencia', s0=20*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s25e1).time_min,
               P=Hydrolysis_file(Files.Small.s25e1).alphaNH(),
               function='p_Valencia', s0=25*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s30e1).time_min,
               P=Hydrolysis_file(Files.Small.s30e1).alphaNH(),
               function='p_Valencia', s0=30*Parameters.Htot, e0=1)

fit_Qi(t=Hydrolysis_file(Files.Small.s40e1).time_min,
               P=Hydrolysis_file(Files.Small.s40e1).alphaNH(),
               function='p_Valencia', s0=40*Parameters.Htot, e0=1)
















# plot_several([Files.Small.s5e1, Files.Small.s10e1,
#             Files.Small.s15e1, Files.Small.s30e1,
#             Files.Small.s40e1], 
#             title='Same enzyme concentration small reactor',
#             legend=['5', '10', '15', '30', '40'])

# plot_several([Files.Large.s10e3, Files.Large.s30e3, 
#             Files.Large.s40e3],
#             title='Same enzyme concentration (3) large reactor',
#             legend=['10', '30', '40'])

# plot_several([Files.Large.s10e1, Files.Large.s10e2, 
#             Files.Large.s10e3],
#             title='Same substrate concentration (10gL) large reactor',
#             legend=['1', '2', '3'])



# plot_several([Files.Large.s10e1, Files.Large.s10e2, Files.Large.s10e3, 
#             Files.Large.s30e3, Files.Large.s40e3],
#             title='Several substrate and enzyme concentrations',
#             legend=['s10e1','s10e2','s10e3', 's30e3', 's40e3'])




# plot_several([Files.Large.s40e3, Files.Large.s40e4, Files.Large.s40e8],
#             title='Same substrate concentration (40gL) large reactor',
#             legend=['3', '4', '8'])


# plot_several([Files.Large.s10e1, Files.Small.s10e1],
#             title='Comparing reactors 10gL - 1',
#             legend=['Large', 'Small'])








##### Uncomment to process the files #####

# Hydrolysis_file(Files.Large.s10e1).save_file_to_model()
# Hydrolysis_file(Files.Large.s10e2).save_file_to_model()
# Hydrolysis_file(Files.Large.s10e3).save_file_to_model()
# Hydrolysis_file(Files.Large.s30e3).save_file_to_model()
# Hydrolysis_file(Files.Large.s40e3).save_file_to_model()