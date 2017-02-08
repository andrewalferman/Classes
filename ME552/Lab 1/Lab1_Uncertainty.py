#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:10:13 2017

This program uses the uncertainties package in Python to automatically
calculate and propagate the uncertainties in the equations for ME 552 Lab 1
Winter 2017.  The method used to propagate the uncertainties is the Kline
McClintock method.

@author: andrewalferman
"""

# Import the required packages
import uncertainties as unc
import numpy as np

"""
The following parameters are for experiment 1.  This program is not configured
to run automatically and solve every set of data simultaneously, so the numbers
below must be manually updated for each of the different measurements with
different flow rates.

NOTE: The values used below may not correspond exactly to any one given run.

Also note that all values will be in meters, Pa, seconds, kg, etc.
"""

# Define the universal constants
# Units of R in J/kg*K
R = 287.05

# Add in the inner diameter of the orifice and the piping upstream
# Will be input in inches and output in meters
intom = 0.0254
d0 = unc.ufloat(0.1248*intom, 0.0004472136*intom)
d1 = unc.ufloat(0.3105*intom, 0.0009617692*intom)

# Add in the pressure measurements upstream and downstream of the orifice
# Will pe input in psi and output in Pascals
psitopa = 6894.76
psioffset = 14.7
P1 = unc.ufloat((psioffset+59.37)*psitopa, 0.25*psitopa)
P2 = unc.ufloat((psioffset+59.43)*psitopa, 0.25*psitopa)

# Add in the temperature measurement upstream of the orifice
# Input and output in Kelvin
T = unc.ufloat(19.08+273.15, 2.2)

# Add in the parameters used to calculate the flow rate from the dry test meter
# V input in scf and output in m^3, t is in seconds
cftocm = 0.0283168
V = unc.ufloat(46.01*cftocm,0.01)
t = unc.ufloat(245,1)

# Rotameter measurements
rot = 7.5 * cftocm
Qobs = unc.ufloat(rot, 0.2125*cftocm)

# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
Q = 60*V/t
A0 = 0.25 * np.pi * d0**2
A1 = 0.25 * np.pi * d1**2
deltaP = np.abs(P1 - P2)
rho1 = P1 / (R * T)
E = 1 / ((1 - (A0/A1)**2)**0.5)
CY = Q / (E * A0 * ((2.*deltaP)/rho1)**0.5)
cfunc = Q.n / Qobs.n
calibratedrot = (cfunc * Qobs) / cftocm

print('Rotameter Observed: {}'.format(Qobs))
print('Discharge coefficient times compressiblity constant: {}'.format(CY))
print('Q Actual: {}'.format(Q))
print('Rotameter correction factor: {}'.format(cfunc))
print('Corrected Rotameter: {}'.format(calibratedrot))