import numpy as np
from process_file import Hydrolysis_file
import matplotlib.pyplot as plt
from parameters import Parameters

def coeff_determination(y_actual, y_pred):
    ymean = np.mean(y_actual)
    sstot = np.sum((y_actual - ymean)**2)
    ssres = np.sum((y_actual - y_pred)**2)

    R2 = 1 - ssres / sstot

    return R2


def X2(y_actual, y_pred, n):
    N = len(y_actual) - n
    ssres = np.sum((y_actual - y_pred)**2)

    return ssres / N

def plot_several(listFiles, title, legend):
    for i in listFiles:
        plt.plot(Hydrolysis_file(i).time_min,
         Hydrolysis_file(i).alphaNH(), 'o', fillstyle='none')
    
    plt.xlabel('time / min')
    plt.ylabel('P / mM')
    plt.title(title)
    plt.legend(legend)
    plt.show()


def a_function(x, k2, Ks, p):

    return k2 * x[1] * x[0] * Parameters.Kp * Ks / \
        (Parameters.Km * Ks * p + Ks * Parameters.Kp * x[0] + Parameters.Kp * x[0]**2)


def b_function(s0, k3):

    Km = Parameters.Km
    Ks = Parameters.Ks
    Kp = Parameters.Kp
    k2 = Parameters.k2
    p = Parameters.p

    num = k3 * Km * Ks * Kp
    den = k2 * (Km * Ks * p + Ks * Kp * s0 + Kp * s0**2)

    return num / den

def plot_final(file, s0, e0):

    time = np.arange(0, 60, 0.1)

    plt.plot(Hydrolysis_file(file).time_min,
         Hydrolysis_file(file).alphaNH(), 'o', fillstyle='none')

    a = a_function([s0, e0], Parameters.k2, Parameters.Ks, Parameters.p)
    b = b_function(s0, Parameters.k3)

    P = (1 / b) * np.log(a * b * time + 1)

    plt.plot(time, P)

    plt.show()
