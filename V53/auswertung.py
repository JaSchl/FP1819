import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.stats import stats

x, y = np.genfromtxt('mode1.txt' , unpack=True)
c, h = np.genfromtxt('mode2.txt' , unpack=True)
j, s = np.genfromtxt('mode3.txt' , unpack=True)

def f(x, m, n, k):
    return m * x**2 + n*x +k
def g(c, a, b, o):
    return a*c**2 + b*c + o
def i(j, d, e, t):
    return d*j**2 + e*j +t
paramsI, covarianceI = curve_fit(f, x, y)
errorsI = np.sqrt(np.diag(covarianceI))
paramsII, covarianceII = curve_fit(g, c, h)
errorsII = np.sqrt(np.diag(covarianceII))
paramsIII, covarianceIII = curve_fit(i, j, s)
errorsIII = np.sqrt(np.diag(covarianceIII))

n = ufloat(paramsI[1], errorsI[1])
m = ufloat(paramsI[0], errorsI[0])
k = ufloat(paramsI[2], errorsI[2])
a = ufloat(paramsII[1], errorsII[1])
b = ufloat(paramsII[0], errorsII[0])
o = ufloat(paramsII[2], errorsII[2])
d = ufloat(paramsIII[1], errorsIII[1])
e = ufloat(paramsIII[0], errorsIII[0])
t = ufloat(paramsIII[2], errorsIII[2])
print(m, n, k)
print(a, b, o)
print(d, e, t)

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
q, l = np.genfromtxt('dämpfungmess.txt' , unpack =True)
r, w = np.genfromtxt('dämpfungtheo.txt' , unpack =True)

plt.figure(2)
plt.plot(q, l, 'r+', label='SWR-Meter')
plt.plot(r, w, 'b+', label='Theoriewerte')
plt.ylabel(r"Dämpfung$\,/\,$dB")
plt.xlabel(r"Mikrometereinstellung$\,/\,$mm")

plt.tight_layout
plt.legend(loc="best")

plt.savefig('dämpfung.pdf')
