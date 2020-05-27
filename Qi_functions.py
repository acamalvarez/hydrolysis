import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy import optimize

from bounds import Bounds
from parameters import Parameters
from files import Files

from process_file import Hydrolysis_file

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



def fit_Qi(t, P, function, s0, e0, symbol='o'):
    '''
    Function for the fitting of the experimental data to 
    the equation of the mechanism from Qi
    '''

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    plt.plot(t, P, symbol, fillstyle='none', markevery=0.05)

    if function == 'ab_Valencia':
        popt_Qi_ab, pcov = optimize.curve_fit(Qi(s0, e0, t).ab_Valencia, t, P, p0=[2, 2],
                                              bounds=(Bounds.lb_ab, Bounds.ub_ab))
        
        plt.plot(t, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab))

        R2 = coeff_determination(P, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab))
        chi = X2(P, Qi(s0, e0, t).ab_Valencia(t, *popt_Qi_ab), 2)

        values = np.zeros(4)

        values[0:2] = popt_Qi_ab
        values[-2] = R2
        values[-1] = chi

        results = {"Parameters":['a', 'b', 'R2', 'X2'], "Values":np.round(values, 3)}

    elif function == 'p_Valencia':

        popt_Qi_p, pcov = optimize.curve_fit(Qi(s0, e0, t).p_Valencia, t, P, p0=[1, 1, 1], bounds=(Bounds.lb_p, Bounds.ub_p))

        plt.plot(t, Qi(s0, e0, t).p_Valencia(t, *popt_Qi_p), 'C2')

        R2 = coeff_determination(P, Qi(s0, e0, t).p_Valencia(t, *popt_Qi_p))
        chi = X2(P, Qi(s0, e0, t).p_Valencia(t, *popt_Qi_p), 3)

        values = np.zeros(5)
        values[0:3] = popt_Qi_p
        values[-2] = R2
        values[-1] = chi

        results = {"Parameters":['k2', 'p', 'k3', 'R2', 'X2'], "Values":np.round(values, 3)}


    elif function == 'sp_Valencia':

        popt_Qi_sp, pcov = optimize.curve_fit(Qi(s0, e0, t).sp_Valencia, t, P, p0=[1, 1, 1, 1], bounds=(Bounds.lb_sp, Bounds.ub_sp))

        plt.plot(t, Qi(s0, e0, t).sp_Valencia(t, *popt_Qi_sp), 'C2')

        print(popt_Qi_sp)

        results = {"Parameters":['k2', 'Ks', 'p', 'k3'], "Values":np.round(popt_Qi_sp, 2)}


    plt.xlabel('Time')
    plt.ylabel('alhpa-NH')
    # plt.legend()
    # plt.show()

    df = pd.DataFrame(results)

    print(df)

    return results["Values"]


def fit_several(files, function, s0s, e0s, reactor):

    symbols = ['o', 's', '^', 'v', '*']

    for i in np.arange(len(files)):
        fit_Qi(t=Hydrolysis_file(files[i]).time_min,
        P=Hydrolysis_file(files[i]).alphaNH(),
        function=function, s0=s0s[i], e0=e0s[i], symbol=symbols[i])


    if reactor == 'large':
    
        plt.legend([r"s$_0$ = 91 mM, e$_0$ = 1 g/L",  '_nolegend_', r"s$_0$ = 91 mM, e$_0$ = 2 g/L", '_nolegend_', 
                r"s$_0$ = 91 mM, e$_0$ = 3 g/L", '_nolegend_',  r"s$_0$ = 273 mM, e$_0$ = 3 g/L", '_nolegend_', 
                r"s$_0$ = 364 mM, e$_0$ = 3 g/L"], loc=2)

    elif reactor == 'small':

        plt.legend([r"s$_0$ = 45 mM",  '_nolegend_', r"s$_0$ = 91 mM", '_nolegend_', 
                r"s$_0$ = 136.5 mM", '_nolegend_',  r"s$_0$ = 273.0 mM"], loc=2)

    plt.xlabel('Time (min)', fontweight='bold', fontsize=12)
    plt.ylabel('P (mM)', fontweight='bold', fontsize=12)
    plt.ylim([0, 16])
    plt.tight_layout()
    plt.savefig('Hydrolysis_' + reactor + '_reactor.svg')