import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
from scipy.stats import stats
import scipy.constants as const
import math

z = 2 *const.e
Ao =16
An = 14
A = 30/2
errorA = 1/2 * np.sqrt((A - Ao)**2+(A-An)**2)
A = unp.uarray(A, errorA)
rho = 1.293
N = rho/ (A* const.u)
Zo = 8
Zn = 7
Z = 15/2
errorZ = 1/2 * np.sqrt((Z - Zo)**2+(Z-Zn)**2)
Z = unp.uarray(Z, errorZ)
v = unp.uarray(1.548, 0.009)*10**7
Io = 13.62
In = 14.53
I = (Io+In)/2
errorI = 1/2 * np.sqrt((I -Io)**2+(I-In)**2)
I = unp.uarray(I, errorI)
print('Mittel I', I)
I =I*const.e

dedx = (4*np.pi*const.e**4 * z**2 * N * Z)/ (const.m_e * v**2 * (4*np.pi*const.epsilon_0)**2) + unp.log((2*const.m_e*v**2)/I)

print('denachdx : ', dedx)
print('Mittel I', I)
print('Mittel A ', A)
print('Mittel Z ', Z)

#----------------------------------------------------------

R = const.R
T = 293.15
M = 28.96 * 10**(-3)  #g/mol

p = (rho * R * T)/M

print('Druck p = ', p)
