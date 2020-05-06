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
        self.e0 = e0
        self.s0 = s0
        self.t_data = t_data


    def ab_Valencia(self, t, a, b):

        return (1 / b) * np.log(a * b * t + 1)



    def a_b(self, t, a, b):

        def dpdt(t, h):

            return a * np.exp(-b * h)

        h = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0],
                      method='RK45', t_eval=self.t_data)

        return h.y[0]

    def no(self, t, k2, k3, Km):
        ''' mechanism with no inhibition
        '''

        def dhdt(t, h):

            a = k2 * self.e0 / self.s0
            b = k3 * Km / k2

            return a * np.exp(-b * h)

        h = solve_ivp(dhdt, [self.t_data[0], self.t_data[-1]], [0],
                      method='RK45', t_eval=self.t_data)

        return h.y[0]

    def s(self, t, k2, k3, Km, Ks):
        ''' mechanism with inhibition by substrate
        '''

        def dhdt(t, h):

            a = k2 * Ks * self.e0 / (self.s0 * Ks + self.s0**2)
            b = k3 * Km * Ks / (k2 * (Ks + self.s0))

            return a * np.exp(-b * h)

        h = solve_ivp(dhdt, [self.t_data[0], self.t_data[-1]],
                      [0], method='RK45', t_eval=self.t_data)

        return h.y[0]

    def p(self, t, k2, k3, Km, Kp, p):
        ''' mechanism with inhibition by product
        '''

        def dhdt(t, h):

            a = k2 * self.e0 * Kp / (self.s0 * Kp + p * Km)
            b = k3 * Km * Kp * self.s0 / (k2 * (self.s0 * Kp + p * Km))

            return a * np.exp(-b * h)

        h = solve_ivp(dhdt, [self.t_data[0], self.t_data[-1]],
                      [0], method='RK45', t_eval=self.t_data)

        return h.y[0]

    def sp(self, t, k2, k3, Km, Kp, p, Ks):
        ''' mechanism with inhibition by substrate and product
        '''

        def dhdt(t, h):

            a = self.s0 * k2 * self.e0 * Ks * Kp / \
                (Ks * Kp * self.s0 + Kp * self.s0**2 + Km * Ks * p)
            b = k3 * Km * Ks * Kp * self.s0 / \
                (k2 * (Ks * Kp * self.s0 + Kp * self.s0**2 + Km * Ks * p))

            return a * np.exp(-b * h)

        h = solve_ivp(dhdt, [self.t_data[0], self.t_data[-1]],
                      [0], method='RK45', t_eval=self.t_data)

        return h.y[0]


def fit_Qi(t, h, function, s0, e0):
    '''
    Function for the fitting of the experimental data to 
    the equation of the mechanism from Qi
    '''

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t, h, 'o', label='Experimental data', fillstyle='none', markevery=0.05)

    if function == 'a_b':
        popt_Qi_ab, pcov = optimize.curve_fit(Qi(e0, s0, t).a_b, t, h, p0=[2, 2],
                                              bounds=(Bounds.lb_ab, Bounds.ub_ab))
        
        ax.plot(t, Qi(e0, s0, t).a_b(t, *popt_Qi_ab),
                'C1-', label='Qi just a and b')

        R2 = coeff_determination(h, Qi(e0, s0, t).a_b(t, *popt_Qi_ab))
        chi = X2(h, Qi(e0, s0, t).a_b(t, *popt_Qi_ab), 2)

        values = np.zeros(4)

        values[0:2] = popt_Qi_ab
        values[-2] = R2
        values[-1] = chi

        results = {"Parameters":['a', 'b', 'R2', 'X2'], "Values":np.round(values, 3)}

    elif function == 'ab_Valencia':
        popt_Qi_ab, pcov = optimize.curve_fit(Qi(e0, s0, t).ab_Valencia, t, h, p0=[2, 2],
                                              bounds=(Bounds.lb_ab, Bounds.ub_ab))
        
        ax.plot(t, Qi(e0, s0, t).ab_Valencia(t, *popt_Qi_ab),
                'C1-', label='Qi just a and b')

        R2 = coeff_determination(h, Qi(e0, s0, t).ab_Valencia(t, *popt_Qi_ab))
        chi = X2(h, Qi(e0, s0, t).ab_Valencia(t, *popt_Qi_ab), 2)

        values = np.zeros(4)

        values[0:2] = popt_Qi_ab
        values[-2] = R2
        values[-1] = chi

        results = {"Parameters":['a', 'b', 'R2', 'X2'], "Values":np.round(values, 3)}


    elif function == 'no':
        popt_Qi_no, pcov = optimize.curve_fit(Qi(e0, s0, t).no, t, h, p0=[1, 1, 43.92],
                                              bounds=(Bounds.lb_no, Bounds.ub_no))
        
        ax.plot(t, Qi(e0, s0, t).no(t, *popt_Qi_no),
                'C1-', label='Qi No inhibition')

        results = {"Parameters":['k2', 'k3', 'Km'], "Values":np.round(popt_Qi_no, 3)}

    elif function == 's':
        popt_Qi_s, pcov = optimize.curve_fit(Qi(e0, s0, t).s, t, h, p0=[1, 1, 43.92, 1],
                                             bounds=(Bounds.lb_s, Bounds.ub_s))

        ax.plot(t, Qi(e0, s0, t).s(t, *popt_Qi_s),
                'C2-', label='Qi substrate inhibition')

        results = {"Parameters":['k2', 'k3', 'Km', 'Ks'], "Values":np.round(popt_Qi_s, 3)}


    elif function == 'p':
        popt_Qi_p, pcov = optimize.curve_fit(Qi(e0, s0, t).p, t, h, p0=[1, 1, 117.45, 22.394, 1],
                                             bounds=(Bounds.lb_p, Bounds.ub_p))

        ax.plot(t, Qi(e0, s0, t).p(t, *popt_Qi_p),
                'C3-', label='Qi product inhibition')

        results = {"Parameters":['k2', 'k3', 'Km', 'Kp', 'p'], "Values":np.round(popt_Qi_p, 3)}


    elif function == 'sp':
        popt_Qi_sp, pcov = optimize.curve_fit(Qi(e0, s0, t).sp, t, h, p0=[1, 1, 43.92, 1, 1, 1],
                                              bounds=(Bounds.lb_sp, Bounds.ub_sp))

        ax.plot(t, Qi(e0, s0, t).sp(t, *popt_Qi_sp), 'C4-',
                label='Qi substrate/product inhibition')

        results = {"Parameters":['k2', 'k3', 'Km', 'Kp', 'p', 'Ks'], "Values":np.round(popt_Qi_sp, 3)}


    ax.set_xlabel('Time')
    ax.set_ylabel('alhpa-NH')
    ax.legend()
    plt.show()

    df = pd.DataFrame(results)

    print(df)

    return results["Values"]
