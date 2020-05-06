import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy import optimize

from bounds import Bounds
from parameters import Parameters
from files import Files

class Function:
    '''
    This class defines all the functions from the model
    of we developed
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

    def Camilo(self, t, kd, k2):

        def dpdt(t, p):
            
            Km = 117.45
            Kp = 22.394

            return (k2 * self.s0 * self.e0 / Km) * np.exp(-kd * p / k2)

        p = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0],
                        method='RK45', t_eval=self.t_data)

        return p.y[0]

    def Tati(self, t, kd, k2):

        def dpdt(t, p):

            Km = 117.45
            Kp = 22.394

            return (k2 * self.e0 / Km) * np.exp(-kd * Km * p**2 / (2 * k2 * Kp))

        p = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0], method='RK45', t_eval=self.t_data)

        return p.y[0]

    def Tati2(self, t, kd, k2):

        def dpdt(t, p):

            Km = 117.45
            Kp = 22.394

            return (k2 * self.e0 * (self.s0) / Km) * np.exp(-kd * p / k2)

        p = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0], method='RK45', t_eval=self.t_data)

        return p.y[0]

    def Tati3(self, t, kd, k2):

        def dpdt(t, p):

            Km = 117.45
            Kp = 22.394

            # return k2 * (self.e0 - (1 - (self.s0 - p) * p / Kp) * self.e0 * np.exp(-kd * p / k2))
            return k2 * p * np.exp(-kd * p)

        p = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0], method='RK45', t_eval=self.t_data)

        return p.y[0]

    def Qi_p(self, t, k2, k3, h):

        def dpdt(t, p):

            Km = 117.45
            Kp = 22.394

            a = k2 * self.e0 * Kp / (self.s0 * Kp + h * Km)
            b = k3 * Km * Kp * self.s0 / (k2 * (self.s0 * Kp + h * Km))
            
            return a * np.exp(-b * p)

        p = solve_ivp(dpdt, [self.t_data[0], self.t_data[-1]], [0], method='RK45', t_eval=self.t_data)

        return p.y[0]




def fit_function(t, h, function, s0, e0):
    '''
    Function for the fitting of the experimental data to 
    the equation of the mechanism from Qi
    '''

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(t, h, 'o', label='Experimental data', fillstyle='none')

    if function == 'Camilo':
        popt_Cam, pcov = optimize.curve_fit(Function(e0, s0, t).Camilo, t, h, p0=[1, 1],
                                              bounds=(Bounds.lb_Cam, Bounds.ub_Cam))

        ax.plot(t, Function(e0, s0, t).Camilo(t, *popt_Cam), 'C4-',
                label='Camilo product inhibition')

        result = popt_Cam

    elif function == 'Tati':
        popt_Tati, pcov = optimize.curve_fit(Function(e0, s0, t).Tati, t, h, p0=[1, 1],
                                              bounds=(Bounds.lb_Tati, Bounds.ub_Tati))

        ax.plot(t, Function(e0, s0, t).Tati(t, *popt_Tati), 'C4-',
                label='Tati product inhibition')

        result = popt_Tati

    elif function == 'Tati2':
        popt_Tati2, pcov = optimize.curve_fit(Function(e0, s0, t).Tati2, t, h, p0=[1, 1],
                                              bounds=(Bounds.lb_Tati, Bounds.ub_Tati))

        ax.plot(t, Function(e0, s0, t).Tati2(t, *popt_Tati2), 'C4-',
                label='Tati product inhibition')

        result = popt_Tati2

        print('kd, k2', result)

    elif function == 'Tati3':
        popt_Tati3, pcov = optimize.curve_fit(Function(e0, s0, t).Tati3, t, h, p0=[1, 1],
                                              bounds=(Bounds.lb_Tati, Bounds.ub_Tati))

        ax.plot(t, Function(e0, s0, t).Tati3(t, *popt_Tati3), 'C4-',
                label='Tati product inhibition')

        result = popt_Tati3

        print('kd, k2', result)

    ax.set_xlabel('Time')
    ax.set_ylabel('alpha-NH')
    ax.legend()
    plt.show()

    return result