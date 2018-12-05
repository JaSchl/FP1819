import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.stats import stats

x, y = np.genfromtxt('mode1.txt' , unpack=True)
c, h = np.genfromtxt('mode2.txt' , unpack=True)
j, s = np.genfromtxt('mode3.txt' , unpack=True)

def f(x, m, n):
    return m * x**2 + n
def g(c, a, b):
    return a*c**2 + b
def i(j, d, e):
    return d*j**2 + e
paramsI, covarianceI = curve_fit(f, x, y)
errorsI = np.sqrt(np.diag(covarianceI))
paramsII, covarianceII = curve_fit(g, c, h)
errorsII = np.sqrt(np.diag(covarianceII))
paramsIII, covarianceIII = curve_fit(i, j, s)
errorsIII = np.sqrt(np.diag(covarianceIII))

n = ufloat(paramsI[1], errorsI[1])
m = ufloat(paramsI[0], errorsI[0])
a = ufloat(paramsII[1], errorsII[1])
b = ufloat(paramsII[0], errorsII[0])
d = ufloat(paramsIII[1], errorsIII[1])
e = ufloat(paramsIII[0], errorsIII[0])
print(m, n)
print(a, b)
print(d, e)

L1plot = np.linspace(50, 70)
L2plot = np.linspace(90, 110)
L3plot = np.linspace(140, 170)

plt.figure(1)
plt.plot(x, y,'r+', label="1. Mode")
plt.plot(c, h,'b+', label="2. Mode")
plt.plot(j, s,'g+', label="3. Mode")
plt.ylabel(r"$A\,/\,V$")
plt.xlabel(r"$U\,/\,V$")
plt.plot(L1plot, f(L1plot, *paramsI) , 'r-')
plt.plot(L2plot, f(L2plot, *paramsII) , 'b-')
plt.plot(L3plot, f(L3plot, *paramsIII) , 'g-')

plt.tight_layout
plt.legend(loc="best")

plt.savefig('moden.pdf')

##################################Dämpfung#############################
#k, l = np.genfromtxt('dämpfung.txt')
#
#plt.plot(k, l, 'r+', label='SWR-Meter')
#plt.ylabel(r"Dämpfung\,/\,dB")
#plt.xlabel(r"Mikrometereinstellung\,/\,mm")
#
#plt.tight_layout
#plt.legend(doc="best")
#
#plt.savefig('dämpfung.pdf')
