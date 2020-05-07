import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy import optimize

from bounds import Bounds
from parameters import Parameters
from files import Files

from functions import coeff_determination, X2

class Qi:
    '''
    This class defines all the functions from the model
    of Qi
    '''

    def __init__(self, s0, e0, t_data):
        '''
        Initializes the functions and the data.
        You must enter the following variables:
        - e0: initial enzyme concentration
        - s0: initial substrate concentration
        - t_data: time data as a numpy array
        '''
        self.s0 = s0
        self.e0 = e0
        self.t_data = t_data


    def ab_Valencia(self, t, a, b):

        return (1 / b) * np.log(a * b * t + 1)

    def p_Valencia(self, t, k2, p, k3):

        a = k2 * self.e0 * self.s0 * Parameters.Kp / (Parameters.Km * p + Parameters.Kp * self.s0)
        b = k3 * Parameters.Km * Parameters.Kp / (k2 * (Parameters.Km * p + Parameters.Kp * self.s0))

        return (1 / b) * np.log(a * b * t + 1)

    def sp_Valencia(self, t, k2, Ks, p, k3):

        a = k2 * self.e0 * self.s0 * Parameters.Kp * Ks / (Parameters.Km * Ks * p + Ks * Parameters.Kp * self.s0 + Parameters.Kp * self.s0**2)
        b = k3 * Parameters.Km * Ks * Parameters.Kp / (k2 * (Parameters.Km * Ks * p + Ks * Parameters.Kp * self.s0 + Parameters.Kp * self.s0**2))

        return (1 / b) * np.log(a * b * t + 1)



def fit_Qi(t, P, function, s0, e0):
    '''
    Function for the fitting of the experimental data to 
    the equation of the mechanism from Qi
    '''

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t, P, 'o', label='Experimental data', fillstyle='none', markevery=0.05)

    if function == 'ab_Valencia':
        popt_Qi_ab, pcov = optimize.curve_fit(Qi(s0, e0, t).ab_Valencia, t, P, p0=[2, 2],
                                              bounds=(Bounds.lb_ab, Bounds.ub_ab))
        
        ax.plot(t, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab),
                'C1-', label='Qi just a and b')

        R2 = coeff_determination(P, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab))
        chi = X2(P, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab), 2)

        values = np.zeros(4)

        values[0:2] = popt_Qi_ab
        values[-2] = R2
        values[-1] = chi

        results = {"Parameters":['a', 'b', 'R2', 'X2'], "Values":np.round(values, 3)}

    elif function == 'p_Valencia':

        popt_Qi_p, pcov = optimize.curve_fit(Qi(s0, e0, t).p_Valencia, t, P, p0=[4.668, 0.786, 2.754], bounds=(Bounds.lb_p, Bounds.ub_p))

        ax.plot(t, Qi(s0, e0, t).p_Valencia(t, *popt_Qi_p), 'C2', label='Valencia product')

        results = {"Parameters":['k2', 'p', 'k3'], "Values":np.round(popt_Qi_p, 2)}


    elif function == 'sp_Valencia':

        popt_Qi_sp, pcov = optimize.curve_fit(Qi(s0, e0, t).sp_Valencia, t, P, p0=[1, 1, 1, 1], bounds=(Bounds.lb_sp, Bounds.ub_sp))

        ax.plot(t, Qi(s0, e0, t).sp_Valencia(t, *popt_Qi_sp), 'C2', label='Valencia substrate and product')

        print(popt_Qi_sp)

        results = {"Parameters":['k2', 'Ks', 'p', 'k3'], "Values":np.round(popt_Qi_sp, 2)}


    ax.set_xlabel('Time')
    ax.set_ylabel('alhpa-NH')
    ax.legend()
    plt.show()

    df = pd.DataFrame(results)

    print(df)

    return results["Values"]
