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



#Foliendicke

p, oF= np.genfromtxt('foliendicke.txt' , unpack=True)
p2, mF= np.genfromtxt('foliendicke2.txt' , unpack=True)

def f(p2, a, b):
        return a*p2 + b
paramsI, covarianceI = curve_fit(f, p2, mF)
errorsI = np.sqrt(np.diag(covarianceI))
a = ufloat(paramsI[0], errorsI[0])
b = ufloat(paramsI[1], errorsI[1])
print('Parameter')
a = a* 10**3
print('a =', a, 'V/mbar')
print('b =', b, 'V')

def f(p, c, d):
        return c*p + d
paramsII, covarianceII = curve_fit(f, p, oF)
errorsII = np.sqrt(np.diag(covarianceII))
c = ufloat(paramsII[0], errorsII[0])
d = ufloat(paramsII[1], errorsII[1])
print('Parameter')
print(c, 'V/mbar')
print(d, 'V')

L1plot = np.linspace(0, 120)
L2plot = np.linspace(0, 150)
plt.figure(1)
plt.plot(p2, mF,'r+', label="Messwerte mit Folie")
plt.plot(L1plot, f(L1plot, *paramsI) , 'm-', label='Regression')
plt.plot(p, oF, 'g+', label='Messwerte ohne Folie')
plt.ylabel(r"$p/ \mathrm{mbar}$")
plt.xlabel(r"$I / \mathrm{V}$")

plt.plot(L2plot, f(L2plot, *paramsII), 'c-', label='Regression')
plt.tight_layout
plt.legend(loc="best")

plt.savefig('plotfoliendicke.pdf')

deltaI = d-b
print('Energieverlust:', deltaI, 'V')

Ealpha = 5.486*10**6 #https://www.phywe.de/de/alpha-und-photodetektor.html
deltaE = (d/deltaI)*Ealpha
print('Differenzenergie:', deltaE, 'eV')

Ealpha = const.e *Ealpha
malpha = 6.64*10**(-27) #https://www.spektrum.de/lexikon/physik/alphastrahlung/393
valpha = np.sqrt((2*Ealpha)/malpha)
#valpha = valpha * ((10**(-3))/(60*60))
print('Geschwindigkeit Alphateilchen:', valpha, 'm/s')
rho = 19.282*(10**(-3)/10**(-6)) #g/cm^3 http://www.periodensystem.info/elemente/gold/
A = 196.97
matom = A *const.u
N = rho / matom
print('N=', N, '1/m^3')

Ienergiegold =9.225 #eV http://www.uniterra.de/rutherford/tab_iong.htm
Ienergiegold = Ienergiegold * const.e
deltaE = deltaE *  const.e
Zhe = 2 #http://www.periodensystem.info/elemente/helium/
Zau = 79 #SIEHE quelle rho, malpha

dx = (deltaE*const.m_e*valpha**2*(4*np.pi*const.epsilon_0)**2)/(4*np.pi*const.e**4*Zhe**2*N*Zau*np.log((2*const.m_e*valpha**2)/(Ienergiegold)))
print('dx=', dx)
#x1 = deltaE*const.e*const.m_e*valpha**2*(4*np.pi*const.epsilon_0)**2
#x2 = 4*np.pi*const.e**4*Zhe**2*N*Zau*np.log(2*const.m_e*valpha**2/(Ienergiegold*const.e))
#x2b = np.log(2*const.m_e*valpha**2/(Ienergiegold*const.e))
#print(x1, x2, x2b)
dxtheo =2*10**(-6)
deltax =(dx-dxtheo)/dxtheo*100
print('prozentuelle Abweichung', deltax, '%')


#Differentieller Streuquerschnitt

theta, c, t = np.genfromtxt('wirkung.txt', unpack=True)
theta = (2*np.pi)/360 *theta
cpros = c/t
errorc =np.sqrt(cpros)
cpros = unp.uarray(cpros, errorc)

x= 2*10**(-3)
y= 10*10**(-3)
L = 4.1*10**(-2)
dOmega = 4* unp.arctan((x)/(2*L)) * unp.arctan(y/(2*L))
print('dOmega =', dOmega)

#zhe= 3.204 * 10**(-19) # =2*e
print('c pro s:', cpros)
#print('Fehler:', errorc)
dsigmadomega2 = 1/((4*np.pi*const.epsilon_0)**2) * ((Zhe* Zau * const.e**2)/(4*Ealpha))**2 * 1/(np.sin(theta/2))**4
dsigmadomega = c/(A * N* dx * dOmega**2)
#print(stds(dsigmadomega))
print('Wirkungsquerschnitt:', dsigmadomega)
print('Wirkungsquerschnitt2:', dsigmadomega2)
#Winkel = np.sin(theta/2)**4
#print('Winkel', Winkel)
plt.figure(2)
plt.errorbar(theta, noms(dsigmadomega), stds(dsigmadomega) , fmt= 'x', label= r'$\frac{d\sigma}{d\Omega}_c$')
plt.plot(theta, dsigmadomega2,'r+', label=r"$\frac{d\sigma}{d\Omega}_s$")
plt.ylabel(r"$ \frac{d\sigma}{d\Omega} / \mathrm{m^{-1}}$")
plt.xlabel(r"$ \theta / \mathrm{rad}$")
plt.grid()
plt.tight_layout
plt.legend(loc="best")
plt.savefig('plotdsigmadomega.pdf')

abweichung = (dsigmadomega -dsigmadomega2)/dsigmadomega2 *100
print('Abweichung Wq', abweichung)

#Mehrfachstreuung

Auvier = 4* 10**(-6)
thetavier =5
cvier =9
tvier= 600
cprosvier = cvier/tvier
print('c pro s 4mu m =', cprosvier )
Auzwei = 2* 10**(-6)
thetazwei = 4
czwei = 1810
tzwei = 300
cproszwei = czwei/tzwei
print('c pro s 2mu m =', cproszwei)  # Mittelwert aller???? => 2,5202 <= bessser!!! eig sowas in Richtung 5,074499

#ZAbhängigkeit

Agold =197
Abis = 209
Aalu= 27

mgold = Agold *const.u
mbis = Abis *const.u
malu = Aalu *const.u

rhogold = 19.282 *10**(-3)/10**(-6)
rhobis = 9.807 *10**(-3)/10**(-6)
rhoalu =  2.698 *10**(-3)/10**(-6)

cgold = unp.uarray(5.074, 0.958)
cbis = unp.uarray(0.35, 0.05)
calu = unp.uarray(0.68, 0.06)

Ngold = rhogold/mgold
Nbis = rhobis/mbis
Nalu = rhoalu/malu
print('Ngba', Ngold, Nbis, Nalu)
#
#dsigmadomegaZgold = cgold/(Agold * Ngold* dx * dOmega**2)
#dsigmadomegaZbis = cbis/(Abis * Nbis* dx * dOmega**2)
#dsigmadomegaZalu = calu/(Aalu * Nalu* dx * dOmega**2)
#
#print('ds dO', dsigmadomegaZgold, dsigmadomegaZbis, dsigmadomegaZalu)
Z, cZ, errorcZ, Ngba= np.genfromtxt('zabhängigkeit.txt' , unpack=True)
cZ = unp.uarray(cZ, errorcZ)
Ngba = Ngba* 10**(28)
dxr = 2*10**(-6)
ZAbhängigkeit= cZ/(Ngba*dxr)
print('krma', ZAbhängigkeit)

plt.figure(3)
plt.errorbar(Z, noms(ZAbhängigkeit), stds(ZAbhängigkeit), fmt= 'x', label= 'Werte')
#plt.plot(Z, dsigmadomega2,'r+', label=r"$\frac{d\sigma}{d\Omega}_s$")
plt.ylabel(r"$ \frac{c}{Ndx} / \mathrm{m^2s^{-1}}$")
plt.xlabel(r"$ Z $")
plt.tight_layout
plt.grid()
plt.legend(loc="best")
plt.savefig('plotzabhängigkeit.pdf')


#Aktivität der Probe
c = 3231
t = 120
A = c/t* (4*np.pi)/dOmega
print('Die Aktivität beträgt:', A, '1/s')
