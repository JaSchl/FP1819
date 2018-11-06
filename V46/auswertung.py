import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.stats import stats
import scipy.constants as const
import math
#a) BFeld

z, B = np.genfromtxt('bfeld.txt' , unpack=True)

def f(z, a, b, c):
    return b* (z-a)**2 + c
paramsI, covarianceI = curve_fit(f, z, B)
errorsI = np.sqrt(np.diag(covarianceI))
a = ufloat(paramsI[0], errorsI[0])
b = ufloat(paramsI[1], errorsI[1])
c = ufloat(paramsI[2], errorsI[2])
print('Parameter')
print(a)
print(b)
print(c)
print('Scheitelpunkt:')
print((a, c))
L1plot = np.linspace(-15, 15)
plt.figure(1)
plt.plot(z, B,'r+', label="Messwerte")
#Achsen beschriften!!!!!!
plt.ylabel(r"$B(z)/ \mathrm{mT}$")
plt.xlabel(r"$z / \mathrm{cm}$")
plt.plot(L1plot, f(L1plot, *paramsI) , 'b-', label='Regression')
plt.tight_layout
plt.legend(loc="best")

plt.savefig('plotbfeld.pdf')

# b) n-dotiertes GaAs
L1, t11, t12 = np.genfromtxt('ndotiertGaAs.txt', unpack=True)

##Drehwinkel
D1= 1.36*10**(-3)
t1= ((t11-t12)/2)*(np.pi/180)
t1= (t1/D1)
print('Drehwinkel:')
print(t1)

# c) hochreines GaAs
L2, t21, t22 = np.genfromtxt('hochreinGaAs.txt', unpack=True)
##Drehwinkel
D2 = 5.11*10**(-3)
t2= ((t21-t22)/2)*((np.pi/180))
t2 = (t2/D2)
print(t2)

#L1= L1**2
#L2 = L2**2
plt.figure(2)
plt.plot(L1, t1, 'r+',label='n-dotiert')
plt.plot(L2, t2, 'b+',label='hochrein')
plt.ylabel(r"${\theta / \mathrm{\frac{1}{m}}}$")
plt.xlabel(r"${\lambda ^2 / \mathrm{\mu m^2}}$")
#plt.plot(L2plot, f(L2plot, *paramsI) , 'b-', label='Regression')
plt.tight_layout
plt.legend(loc="best")

plt.savefig('plotGaAs.pdf')

#Bestimmung der effektiven Masse.
tn = t1-t2
print("thetaDifferenz:")
print(tn)

l, t = np.genfromtxt('tdiff.txt' , unpack=True)
#t = np.pi/180

def f(tn, m, n):
    return m*tn + n
paramsI, covarianceI = curve_fit(f, l, tn)
errorsI = np.sqrt(np.diag(covarianceI))

m = ufloat(paramsI[1], errorsI[1])
n = ufloat(paramsI[0], errorsI[0])
print('Parameter')
print(m)
print(n)
L2plot = np.linspace(0, 6)
plt.figure(3)
plt.plot(l, tn,'r+', label="Messwerte")
#Achsen beschriften!!!!!!
plt.ylabel(r"${\theta_{diff} / \mathrm{\frac{1}{m}}}$")
plt.xlabel(r"${\lambda^2 / \mathrm{\mu m^2}}$")
plt.plot(L2plot, f(L2plot, *paramsI) , 'b-', label='Regression')
plt.tight_layout
plt.legend(loc="best")

plt.savefig('tdiff.pdf')

## effektive Masse
n2 = 3.6
N = 1.2*10**(10)
c = c*10**(-3)

eps= 8.85*10**(-12)
Meff = (1*(const.e**3*N*c)/(8*np.pi**2*eps*const.c**3 *n*n2))**(0.5)
print()
print(const.e, N, c, const.c, n, n2)
print()
print('Meff:')
print(Meff)
