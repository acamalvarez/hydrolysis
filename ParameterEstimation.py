from scipy import optimize
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def equation(S, Vmax, Km):

    return Vmax * (1 / (1 + Km / S))

def calculate_S_V(x, y):

    popt, pcov = optimize.curve_fit(equation, x, y, p0=[2, 100])
    
    result = {"Parameter":['Vmax','Km'], "Value":np.round(popt, 3)}

    df = pd.DataFrame(result)

    print(df)

S_0 = np.array([0, 27.30, 91.0, 136.5, 182.0, 227.5, 364.0])
v_0 = np.array([0, 0.31, 0.99, 1.48, 1.5921, 1.6841, 1.7428])

S_4 = np.array([0, 27.3, 72.8, 91.0, 136.5, 182.0, 354.9])
v_4 = np.array([0, 0.2, 0.67, 0.88, 1.25, 1.5, 1.66])

S_20 = np.array([0, 45.5, 91.0, 182.0, 227.5, 364.0])
v_20 = np.array([0, 0.31, 0.67, 1.42, 1.5, 1.61])

calculate_S_V(S_0, v_0)
calculate_S_V(S_4, v_4)
calculate_S_V(S_20, v_20)