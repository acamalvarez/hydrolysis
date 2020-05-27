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


def plot_several(listFiles, title, legend, save=""):
    for i in listFiles:
        plt.plot(Hydrolysis_file(i).time_min,
         Hydrolysis_file(i).alphaNH(), 'o', fillstyle='none',
         markevery=0.05)
    
    plt.xlabel('Time (min)', fontweight='bold')
    plt.ylabel('P (mM)', fontweight='bold')
    plt.legend(legend, loc=2)
    plt.ylim([0, 37])
    plt.tight_layout()

    if save != "":
        plt.savefig(save)

    else: plt.show()


def a_function(x, k2, Ks, p):

    # return k2 * x[1] * x[0] * Parameters.Kp * Ks / \
    #     (Parameters.Km * Ks * p + Ks * Parameters.Kp * x[0] + Parameters.Kp * x[0]**2)

    return  k2 * x[1] * x[0] * Parameters.Kp / (Parameters.Km * p + Parameters.Kp * x[0])


def b_function(s0, k3):

    # Km = Parameters.Km
    # Ks = Parameters.Ks
    # Kp = Parameters.Kp
    # k2 = Parameters.k2
    # p = Parameters.p

    # num = k3 * Km * Ks * Kp
    # den = k2 * (Km * Ks * p + Ks * Kp * s0 + Kp * s0**2)

    # return num / den

    return k3 * Parameters.Km * Parameters.Kp / (Parameters.k2 * (Parameters.Km * Parameters.p + Parameters.Kp * s0))

def plot_final(files, s0s, e0, legends):

    fig = plt.figure(num='Hydrolysis curves with the same parameters')
    ax = fig.add_subplot(111)

    symbols = ['o', 's', '^', 'v', '*']
    # markers_on = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

    for i in np.arange(len(files)):

        time = Hydrolysis_file(files[i]).time_min

        ax.plot(Hydrolysis_file(files[i]).time_min,
            Hydrolysis_file(files[i]).alphaNH(), symbols[i], fillstyle='none', label=legends[i])

        a = a_function([s0s[i], e0], Parameters.k2, Parameters.Ks, Parameters.p)
        b = b_function(s0s[i], Parameters.k3)

        P = (1 / b) * np.log(a * b * time + 1)

        ax.plot(time, P)

        R2 = coeff_determination(Hydrolysis_file(files[i]).alphaNH(), P)

        print('R^2 =', np.round(R2, 3))

    ax.set_xlabel('Time (min)', fontweight='bold', fontsize=12)
    ax.set_ylabel('P (mM)', fontweight='bold', fontsize=12)
    ax.legend()
    plt.tight_layout()

    plt.savefig('Hydrolysis curves.svg')



