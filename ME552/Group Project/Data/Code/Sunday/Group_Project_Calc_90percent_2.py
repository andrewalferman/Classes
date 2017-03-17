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
This program is not configured to run automatically and solve every set of data
simultaneously, so the numbers below must be manually updated for each of the
different measurements.

Also note that all values will be in meters, Pa, seconds, kg, etc.
"""

# Define the universal constants
R = 8.3145 # J/mol*K
R_air = 287.058 # J/kg*K, specific gas constant of air
NA = const.Avogadro
hv = 2258000. # J/kg # Heat of vaporization, assumed to be constant
cpw = 4186. # J/kg # Specific heat of water at ambient temperature
cpm = 910. # J/kg*K
Uethanol = 26400000. # J/kg # Energy in ethanol
Methanol = 0.04607 # kg/mol # Density of ethanol
ethrho100 = 0.78934 # kg/L, 100 percent ethanol
ethrho90 = 0.81770 # kg/L, 90.1 percent ethanol
ethrho80  = 0.84270 # kg/L, 80.3 percent ethanol
ethrho70 = 0.86742 # kg/L, 70.1 percent ethanol

# Conversion factors for handy dandy use later, if needed
intom = 0.0254 # Inches to meters
psitopa = 6894.76 # psi to Pascals
psioffset = 14.7 # psig to psia
cftocm = 0.0283168 # cubic feet to cubic meters
inwtopa = 248.84 # Inches of water to Pascals
cfmtocms = 0.000471947 # Cubic feet per minute to cubic meters per second
mbartopa = 100. # Millibars to Pascals
mmtom = .001 # Millimenters to meters
degftometric = 9/5 # Convert F to C or K, not including absolute scale

# Put all the equipment tolerances in here because it's the cool thing to do
thermocouples = 2.2 # Degrees C, but could be 0.75% reading
thermometer = 0.2 * degftometric # Degrees F
barometer = 1. # Millibar
caliper = 0.5 # Millimeters
scale1 = 0.001 # 1 gram accuracy
press = 0.001 # in of H2O
stopwatch = 0.01 # Photograph the stopwatch to achieve this
CO = 0.5 # parts per million
CO2 = 50 # parts per million, may be 5% of reading though
graduatedcyl = 0.0005 # Half a mililiter

"""
MEASUREMENT VALUES START HERE -------------------------------------------------
"""

# Add in all of the measured masses
# All units will be in kg
fuelmass = 0.041 # Measured value
potmeasurement = 0.099 # Measured value
watermassmeasurement = 0.590 # Measured value
ethanolmass = 0.274 # Measured value
bottletotal = 0.304 # Measured value
hotwatermass = 0.435 # Measured value

# Poured volume of fuel in liters
fuelvolume = 0.050 # Measured value

# Temperature measurements
starttemp = 296.25 # Measured value
ambienttemp = 297.26 # Measured value
finaltemp = 373.15 # Solution was boiling, more accurate than thermcouples

# Add in the pressure measurements of the sampling tube
# Will be input and output in inches of water
pressurediff = 0.455 # Measured value

# Absolute pressure in Pascals
barpressure = 1024. # Measured value

# Exhaust gas constituent measurements
COlevel = 15.93 # Measured value, points 120 to 600
CO2level = 295.5 # Measured value, points 120 to 600

# Add in burn time, in seconds
minutes = 24. # Measured value
seconds = 14. # Measured value

"""
CALCULATION STARTS HERE -------------------------------------------------------
"""

# Create ufloats for all of the measured values
fuelv = unc.ufloat(fuelvolume, graduatedcyl)
methanol = unc.ufloat(ethanolmass, scale1)
mbottle = unc.ufloat(bottletotal, scale1)
mpot0 = unc.ufloat(potmeasurement, scale1)
mpot1 = unc.ufloat(watermassmeasurement, scale1)
mpot2 = unc.ufloat(hotwatermass, scale1)
deltaP = unc.ufloat(pressurediff, press)
Pabsolute = unc.ufloat(barpressure, barometer) * mbartopa
COppm = unc.ufloat(COlevel, CO)
CO2ppm = unc.ufloat(CO2level, CO2)
burntime = unc.ufloat(minutes*60. + seconds, 1)
T0 = unc.ufloat(starttemp, thermocouples)
Tambient = unc.ufloat(ambienttemp, thermometer)
TE = unc.ufloat(finaltemp, thermocouples)

# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
mfuel = fuelv * ethrho90
rho_air = Pabsolute / (R_air * Tambient)
mwater1 = mpot1 - mpot0
mwater2 = mpot2 - mpot0
mevap = mpot1 - mpot2
Cx = methanol / mbottle
Ux = Uethanol * Cx
Ereleased = mfuel * Ux
EH2Oevap = mevap * hv
deltaT = TE - T0
Vdottube = (215.*(deltaP**0.5)) * cfmtocms
C = unc.ufloat(1.5998468653396583, 0.02159726001744998) # Calibrated value
Vdot = Vdottube * C
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
print('Relative uncertainty: +/- {:.2f}%'.format((hc.s / hc.n) *100.))
print('-----------------------------')
print('CO PRODUCTION')
print('Nominal value: {:.2f}%'.format(epsilon.n * 100.))
print('Relative uncertainty: +/- {:.2f}%'.format(
        (epsilon.s / epsilon.n) *100.))
