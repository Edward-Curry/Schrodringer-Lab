# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:19:13 2023

@author: ecurry
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import math

e_n_array = []
n = 1
N = 1000
l = 1/(N-1)
psi = np.zeros(N)
y_sqr = 1000
e = (((n**2)*(np.pi**3))/y_sqr) -1.1


xv = np.arange(0,1,1/N)
v_x = (8*(xv-0.5)**2)-1
k_sqr_arr = y_sqr*(e-v_x)
e_arr = []
psie = []
    




def comp_psi(e):
    psi = np.zeros(N)
    psi[0]=0
    psi[1]=1e-4
    
    k_sqr_arr = y_sqr * (e - v_x)
    
    
    for n in range(1,N-1):
        a = (2*(1-((5/12)*(l**2)*(k_sqr_arr[n])))*psi[n])
        b = ((1+((1/12)*(l**2)*(k_sqr_arr[n-1])))*psi[n-1])
        c = (1+((1/12)*(l**2)*(k_sqr_arr[n+1])))
        psi_n1 = (a - b)/c
        psi[n+1] = psi_n1 #psi_n1(n+1)
    
    return psi


while e < 4.5:
    de = 1e-2
    e = e+0.1
    tol = 10e-6
    
    while abs(de) > tol:
        
        psi1 = comp_psi(e)[-1]
        psi2 = comp_psi(e + de)[-1]
        
        
        if psi1*psi2 > 0:
            de=de
            e = e+de
    
        elif psi1*psi2 < 0:
            e = e+de
            de = -de/2
    
    e_n_array.append(e)
    print(e)
    
lower_limit = 0
upper_limit = 1
l = 1.0/(N-1)


result = integrate.simps(comp_psi(e_n_array[n-1])**2, dx = l)
normalised = (comp_psi(e_n_array[n-1])/math.sqrt(result))
int_normalised = integrate.simps(normalised**2, dx = l)
print(int_normalised)

e = (((n**2)*(np.pi**3))/y_sqr) -1
print(e)

x = np.arange(0,1,1/N)
plt.plot(x,normalised)

dd_normalise = np.zeros(N)

for n in range (1,len(normalised)-1):
    dd_normalise[n] =  (normalised[n-1]-(2*normalised[n]) + normalised[n+1])/(l**2)

    
dd_normalise[-1] = 0

result_p = integrate.simps(-1*normalised*dd_normalise, dx = l)

D_p = math.sqrt(abs(result_p))
print("<p^2>", result_p)
print("DELTA p =", D_p)

result_x = integrate.simps(x*normalised**2, dx = l)
result_xx =  integrate.simps(x**2*normalised**2, dx = l)
D_x = math.sqrt(result_xx - (result_x)**2)
print("DELTA X =",D_x)
print("Delta X * DELTA p =",D_x*D_p)

DXDP = np.zeros(10)

for i in range(1,10):
    lower_limit = 0
    upper_limit = 1
    l = 1.0/(N-1)
    
    
    result = integrate.simps(comp_psi(e_n_array[i-1])**2, dx = l)
    normalised = (comp_psi(e_n_array[i-1])/math.sqrt(result))
    int_normalised = integrate.simps(normalised**2, dx = l)
    print(int_normalised)
    
    e = (((i**2)*(np.pi**3))/y_sqr) -1
    print(e)
    
    x = np.arange(0,1,1/N)
    plt.plot(x,normalised)
    plt.xlabel("x")
    plt.ylabel("F(x)")
    plt.title(f"Normalised Wavefunction Harmonic Potential $\epsilon_{i}$")
    filename = f"normalised_wavefunction_plot_Harmonic{i}.pdf"
    plt.savefig(filename)
    plt.show()
    dd_normalise = np.zeros(N)
    
    for n in range (1,len(normalised)-1):
        dd_normalise[n] =  (normalised[n-1]-(2*normalised[n]) + normalised[n+1])/(l**2)
    
        
    dd_normalise[-1] = 0
    
    result_p = integrate.simps(-1*normalised*dd_normalise, dx = l)
    D_p = math.sqrt(abs(result_p))
    print("<p^2>", result_p)
    print("DELTA p =", D_p)
    
    result_x = integrate.simps(x*normalised**2, dx = l)
    result_xx =  integrate.simps(x**2*normalised**2, dx = l)
    D_x = math.sqrt(result_xx - (result_x)**2)
    print("DELTA X =",D_x)
    print("Delta X * DELTA p =",D_x*D_p)
    DXDP[i] = D_x*D_p
    
plt.plot(DXDP)
plt.xlabel("$\epsilon_n$")
plt.ylabel("$\Delta X \Delta P$")
plt.title("$\epsilon_n$ vs $\Delta X \Delta P$")
plt.savefig("DXDP_harmonci.pdf")
plt.show()
e_n_arr_shift = np.roll(e_n_array, 1)
e_dif = abs(e_n_array-e_n_arr_shift)
e_difference = np.roll(e_dif,-1)
print(e_difference)
n_arr = []
log_arr=[]
for q in range(1,19):
    n_arr.append(math.log(q))
    log_arr.append(math.log(e_difference[q-1]))
    
plt.plot(n_arr,log_arr)
plt.xlabel("log(n)")
plt.ylabel("log(adjavecnt energery level difference")
plt.title("log(n) vs  log(energy level difference)")
plt.savefig("energy level difference.pdf")
plt.show()