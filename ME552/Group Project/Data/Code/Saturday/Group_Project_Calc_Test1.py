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
cpm = 490. # J/kg*K
Uethanol = 26400000. # J/kg # Energy in ethanol
Methanol = 0.04607 # kg/mol # Density of ethanol

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

# Add in all of the measured masses
# All units will be in kg
canmass = 0.011 # Measured value
fuelmass = 0.046 # Measured value
potmeasurement = 0.099 # Measured value
watermassmeasurement = 0.246 # Measured value
ethanolmass = 1. # Dummy values because 100% ethanol was used
bottletotal = 1. # Dummy values because 100% ethanol was used
hotwatermass = 0.122 # Measured value

# Create ufloats for all of the measured values
mcan0 = unc.ufloat(canmass, scale1)
mcan1 = unc.ufloat(fuelmass + canmass, scale1)
methanol = unc.ufloat(ethanolmass, scale1)
mdiluent = unc.ufloat(bottletotal - ethanolmass, scale1)
mpot0 = unc.ufloat(potmeasurement, scale1)
mpot1 = unc.ufloat(watermassmeasurement + potmeasurement, scale1)
mpot2 = unc.ufloat(hotwatermass + potmeasurement, scale1)

# Note the poured volume of fuel
fuelvolume = 0.060 # Measured volume in liters
# Create ufloat for volume, will be turned into mass measurement
fuelv = unc.ufloat(fuelvolume, graduatedcyl)

# Add in the pressure measurements of the sampling tube
# Will be input and output in inches of water
deltaP = unc.ufloat(0.445, press) # NEED TO ADD

# Absolute pressure in Pascals
Pabsolute = unc.ufloat(1018.5, barometer) * mbartopa # Measured value

# Add in the CO and CO2 readings
COppm = unc.ufloat(10., CO) # NEED TO ADD
CO2ppm = unc.ufloat(200., CO2) # NEED TO ADD

# Add in burn time, in seconds
minutes = 24. # NEED TO ADD
seconds = 50. # NEED TO ADD
burntime = unc.ufloat(minutes*60. + seconds, 1) # Measured value

# Add in the temperature measurements
# Input and output in Kelvin
T0 = unc.ufloat(294.65, thermocouples) # NEED TO ADD
Tambient = unc.ufloat(297.65, thermometer) # Measured value
TE = unc.ufloat(372.65, thermocouples) # NEED TO ADD

# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
rho_air = Pabsolute / (R_air * Tambient)
mwater1 = mpot1 - mpot0
mwater2 = mpot2 - mpot0
mevap = mwater1 - mwater2
mfuel = mcan1 - mcan0
Cx = methanol / (methanol + mdiluent)
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
