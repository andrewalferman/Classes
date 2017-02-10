#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 10:10:13 2017

This program uses the uncertainties package in Python to automatically
calculate and propagate the uncertainties in the equations for the ME 552
Winter 2017 Group 1 project.  The method used to propagate the uncertainties is
the Kline McClintock method.

@author: andrewalferman
"""

# Import the required packages
import uncertainties as unc
import scipy.constants as const

"""
The following parameters are for experiment 1.  This program is not configured
to run automatically and solve every set of data simultaneously, so the numbers
below must be manually updated for each of the different measurements.

Also note that all values will be in meters, Pa, seconds, kg, etc.
"""

# Define the universal constants
R = 8.3145 # J/mol*K
NA = const.Avogadro
hv = 2258000. # J/kg
cpw = 4186. # J/kg
cpm = 490. # J/kg*K
Uethanol = 26400000. # J/kg
Methanol = 0.04607 # kg/mol

# Conversion factors for handy dandy use later, if needed
intom = 0.0254 # Inches to meters
psitopa = 6894.76 # psi to Pascals
psioffset = 14.7 # psig to psia
cftocm = 0.0283168 # cubic feet to cubic meters
inwtopa = 248.84 # Inches of water to Pascals
cfmtocms = 0.000471947 # Cubic feet per minute to cubic meters per second

# Put all the equipment tolerances in here because it's the cool thing to do
thermocouples = 2.2 # degrees C, but could be 0.75% reading
scale1 = 0.001 # 1 gram accuracy
scale2 = 0.0001 # If we can get 0.1 gram accuracy out of the better scale
pressscale = 10. # Assume that the range will be 10 inches H2O
press = 0.00006 * pressscale
dtm = 0.01 # Dry test meter, in cfm
stopwatch = 0.01 # Photograph the stopwatch to achieve this
CO = 0.5 # parts per million
CO2 = 50 # parts per million, may be 5% of reading though
dpt = 50 # Pascals, if we even use this
tempsensor = 0.5 # Deg C, if we use this
dpgage = 1.*0.02*inwtopa # If we even use this gage from the PEMS

# Add in all of the measured masses
# All units will be in kg
mcan0 = unc.ufloat(0.010, scale1) # Similar to values found at Aprovecho
mcan1 = unc.ufloat(0.050, scale1) # Similar to values found at Aprovecho
methanol = unc.ufloat(0.95, scale1) # Assuming 95% concentration
mdiluent = unc.ufloat(0.050, scale1) # Assuming 95% concentration
mpot0 = unc.ufloat(0.5, scale1) # Assuming the pot will weigh around 0.5 kg
mpot1 = unc.ufloat(1.0, scale1) # Assuming 0.5 kg / 0.5 L water added
mpot2 = unc.ufloat(0.995, scale1) # Assuming we lose 5g of water due to boiling

# Add in the pressure measurements of the sampling tube
# Will be input and outputin inches of water
Pstag = unc.ufloat(0.1, press) # The one meter went up to 1 in, so 0.1 in
                               # sounds at least somewhat reasonable.
Pstat = unc.ufloat(0, press)    # I think this should be 0 psig.

# Absolute pressure in Pascals
Pabsolute = unc.ufloat(101300, 0.0007 * 101300)

# Add in the CO and CO2 readings
COppm = unc.ufloat(10., CO)
CO2ppm = unc.ufloat(200., CO2)

# Add in burn time, in seconds
burntime = unc.ufloat(15.*60., stopwatch) # Assuming that it will burn 15 min

# Add in the temperature measurements
# Input and output in Kelvin
T0 = unc.ufloat(292.23, thermocouples) # This value was what we had in lab 1
Tambient = T0 # Should equilibriate to this.
TE = unc.ufloat(373.15, thermocouples) # Assuming boiling point at sea level


# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
mwater1 = mpot1 - mpot0
mwater2 = mpot2 - mpot0
mevap = mwater1 - mwater2
mfuel = mcan1 - mcan0
Cx = methanol / (methanol + mdiluent)
Ux = Uethanol * Cx
Ereleased = mfuel * Ux
EH2Oevap = mevap * hv
deltaT = TE - T0
deltaP = Pstag - Pstat
Vdot = 215.*(deltaP**0.5) * cfmtocms
V = Vdot * burntime
COpartial = Pabsolute * COppm / 1.e6
CO2partial = Pabsolute * CO2ppm / 1.e6
nco = (NA * COpartial * V) / (R * Tambient)
nco2 = (NA * CO2partial * V) / (R * Tambient)
nethanol = NA * mfuel * Cx / Methanol
EH2Oheat = mwater2*cpw*deltaT + mpot0*cpm*deltaT
hc = (EH2Oheat + EH2Oevap) / Ereleased
epsilon = nco / (2 * nethanol)

print('-----------------------------')
print('THERMAL EFFICIENCY (h_c)')
print('Nominal value: {:.2f}%'.format(hc.n * 100.))
print('Uncertainty: +/- {:.2f}%'.format(hc.s *100.))
print('-----------------------------')
print('CO PRODUCTION')
print('Nominal value: {:.2f}%'.format(epsilon.n * 100.))
print('Uncertainty: +/- {:.2f}%'.format(epsilon.s *100.))
