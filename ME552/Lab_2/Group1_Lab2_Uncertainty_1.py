#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:47:49 2017

This program uses the uncertainties package in Python to automatically
calculate and propagate the uncertainties in the equations for the ME 552
Winter 2017 Lab 2.  The method used to propagate the uncertainties is
the Kline McClintock method.

@authors: Andrew Alferman and Matt Zaiger
"""

# Import the required packages
import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.pyplot as plt
import scipy.constants as const
import os
import powerlaw as pl

import reader as r

# Define the universal constants
R = 286.9 # J/kg*K

# For the purposes of the equation, we make an assumption about n
n = 1./3.

### Call reader function ###
V, u_stds, Tambreadings, Pressreadings = r.reader()

# List out all of the experimental data needed for the calculations
OHR = 1.22

# Conversion factors for handy dandy use later, if needed
intom = 0.0254 # Inches to meters
psitopa = 6894.76 # psi to Pascals
cftocm = 0.0283168 # cubic feet to cubic meters
inwtopa = 248.84 # Inches of water to Pascals
cfmtocms = 0.000471947 # Cubic feet per minute to cubic meters per second

# Put all the equipment tolerances in here because it's the cool thing to do
thermocouples = 1.0 # degrees C (Type T), but could be 0.75% reading
pressscale = 10.0 # Assume that the range will be 10 inches H2O
press = 0.00060 * pressscale
# hwrp = Hotwire resistance, percent
hwrp = 0.001 # Table states "0.1% +/- .01 Ohms" so 2 numbers are used
# hwra = Hotwire resistance, absolute
hwra = 0.01 # Ohms
hotwireoffset = 0.0015 # Percent accuracy
hotwiregain = 0.0015 # Percent accuracy

uv = ((hotwireoffset*V)**2 + (V*hotwiregain)**2 + u_stds**2)**0.5

Volts = [unc.ufloat(V[i], uv[i]) for i in range(len(V))]

# Add in the pressure measurements
# Input in inches of water, output in Pascals
P_atm = unc.ufloat(100200, 100)
P_diffs = [unc.ufloat(p, press*inwtopa) for p in Pressreadings]

# Add in the temperature measurements
# Input and output in Kelvin
T_atms = [unc.ufloat(t, thermocouples) for t in Tambreadings]

# Add in the resistances
# Input and output in ohms
R_0c = 5.85
R_0c = unc.ufloat(R_0c, R_0c * hwrp + hwra)
R_100c = 1.33 + R_0c
R_int = 0.50
R_int = unc.ufloat(R_int, R_int * hwrp + hwra)
R_op = 7.44
R_op = unc.ufloat(R_op, R_op * hwrp + hwra)

# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
rho_airs = [P_atm / (R * T_atm) for T_atm in T_atms]
upressures = [((2*P_diffs[i])/rho_airs[i])**0.5 for i in range(len(rho_airs))]

T_w = 100 * ((R_op - R_0c) / (R_100c - R_0c))


A = (Volts[0]**2) / (R_op * (T_w - T_atms[0]))

B = [None]*len(Volts)
u = [None]*len(Volts)

for i in range(len(Volts)):
    if i == 0:
        B[i]= unc.ufloat(-0.001, 0.012)
        u[i]= unc.ufloat(0.1, 0.01)
    else:
        B[i] = (((Volts[i]*Volts[i])/(R_op*(T_w-T_atms[i]))) - A) /\
                 (upressures[i]*upressures[i])
        u[i] = ((((Volts[i]*Volts[i])/(R_op*(T_w-T_atms[i])))-A)/B[i])**0.5


totalbs = 0.0
for i in range(len(B)):
    totalbs += B[i].s
Bsavg = totalbs/len(B)


Bnorm = [val.n for val in B]
Bavg = min(Bnorm)
Anorm = A.n

umadeup = np.linspace(0, 10)
#voltage = [val.n for val in V]




############################ Figure Stuff #####################################
#print(u)
#print(V)

savepath = '/nfs/stak/students/a/alfermaa/Desktop/Classes/ME552/Lab_2/Data/'


plt.figure(figsize = (8,6), dpi = 300)
unorm = np.array([val.n for val in u])
#check = pl.Fit(unorm,V)
#### OHR122
a122 = 0.128165
b122 = 1.2324
n122 = 0.6148
v122 = a122 + b122*umadeup**n122

#### OHR139
a139 = 0.0
b139 = 0.1975
n139 = 0.5079
v139 = a139 + b139*umadeup**n139
#### OHR153
a153 = 0.013131
b153 = 1.7954
n153 = 0.4928
v153 = a153 + b153*umadeup**n153
#print check.alpha
us = [val.s for val in u]
plt.errorbar(V, unorm, yerr = us, fmt ='o', label = 'Data')
z = np.polyfit(V, unorm, 2)
p = np.poly1d(z)
plt.plot(V, p(V), 'r--', label = 'Polynomial Fit')
plt.plot(v122, umadeup,'g-.', label = '$v = A + B*u^n$')
plt.title('Calibration Curve [OHR = 1.22]',fontsize=14)
plt.xlabel('Voltage (V)', fontsize=12)
plt.ylabel('Velocity (m/s)', fontsize=12)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(loc = 2, fontsize = 14)

#plt.plot(a[0], 'ro' )
#plt.show()    #graph()
plt.savefig(savepath + 'OHR122_cal')
#print("y=(%.6f)x^2+(%.6f)x+(%.6f)"%(z[0],z[1], z[2]))

a122 = unc.ufloat(a122, A.s)
b122 = unc.ufloat(b122, Bsavg)

plt.figure(figsize = (8,6), dpi = 300)
Sen122 = (1/n122) * ((v122 - a122)/b122)**((1/n122) - 1) * (1/b122)
Sen139 = (1/n139) * ((v139 - a139)/b139)**((1/n139) - 1) * (1/b139)
Sen153 = (1/n153) * ((v153 - a153)/b153)**((1/n153) - 1) * (1/b153)
v122err = v122 * hotwiregain
v139err = v139 * hotwiregain
v153err = v153 * hotwiregain
Sen122err = unp.std_devs(Sen122)
Sen139err = unp.std_devs(Sen139)
Sen153err = unp.std_devs(Sen153)
Sen122 = unp.nominal_values(Sen122)
Sen139 = unp.nominal_values(Sen139)
Sen153 = unp.nominal_values(Sen153)
plt.errorbar(v122, Sen122, yerr = Sen122err, xerr= v122err, fmt='o', label='OHR = 1.22')
plt.errorbar(v139, Sen139, yerr = Sen139err, xerr= v122err, fmt='o', label='OHR = 1.39')
plt.errorbar(v153, Sen153, yerr = Sen153err, xerr= v122err, fmt='o', label='OHR = 1.53')
plt.xlabel('Voltage (V)', fontsize=12)
plt.ylabel('Sensitivity (m/sV)', fontsize=12)
plt.title('Voltage Sensitivity ($S_u$)',fontsize=14)
plt.legend(fontsize = 14)
#plt.ylim(0,5)
plt.savefig(savepath + 'Volt_Sensitivity')

#u_rms = 0.
#for i in u_stds:
#    u_rms += i
#u_rms = np.sqrt(u_rms)
#Y = u_rms / u.n

#print('-----------------------------')
#print('Flow Velocity')
#print('Nominal value: {:.2f}%'.format(u.n))
#print('Relative uncertainty: +/- {:.2f}%'.format((u.s / u.n) *100.))
#print('-----------------------------')
#print('Turbulence Intensity')
#print('Nominal value: {:.2f}%'.format(Y.n))
#print('Relative uncertainty: +/- {:.2f}%'.format((Y.s / Y.n) *100.))
