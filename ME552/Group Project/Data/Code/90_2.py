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
import scipy.constants as const
import pint

"""
This program is not configured to run automatically and solve every set of data
simultaneously, so the numbers below must be manually updated for each of the
different measurements.

Also note that all values will be in meters, Pa, seconds, kg, etc.
"""

# Define the universal constants
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
R = 8.3145 * ureg.joule / (ureg.mol * ureg.kelvin)
# Specific gas constant of air
R_air = 287.058 * ureg.joule / (ureg.kilogram * ureg.kelvin)
NA = (const.Avogadro / ureg.mol).plus_minus(0)
# Heat of vaporization, assumed to be constant
hv = 2258000. * ureg.joule / ureg.kilogram
# Specific heat of water at ambient temperature
cpw = 4186. * ureg.joule / (ureg.kilogram * ureg.kelvin)
# Specific heat of aluminum at ambient temperature
cpm = 910. * ureg.joule / (ureg.kilogram * ureg.kelvin)
# Energy density of ethanol
Uethanol = 29.64 * ureg.megajoule/ureg.kilogram
# Density of ethanol
Methanol = 0.04607 * ureg.kilogram / ureg.mol
# Density of 100 percent ethanol
ethrho100 = 0.78934 * ureg.kilogram / ureg.liter
# Density of 90.1 percent ethanol
ethrho90 = 0.81770 * ureg.kilogram / ureg.liter
# Density of 80.3 percent ethanol
ethrho80  = 0.84270 * ureg.kilogram / ureg.liter
# Density of 70.1 percent ethanol
ethrho70 = 0.86742 * ureg.kilogram / ureg.liter
# Molecular weight of CO
MCO = 28.01 * ureg.gram / ureg.mol

# Put all the equipment tolerances in here
thermocouples = 2.2 * ureg.kelvin # Based on Type K thermocouple documentation
thermometer = 0.2 * ureg.rankine # Based on the resolution of the thermomenter
barometer = 0.5 * ureg.millibar # Based on the resolution of the barometer
caliper = 0.5 * ureg.millimeters # Based on the resolution of the caliper used
scale1 = 1 * ureg.gram # Based on the scale documentation
press = 0.005 * ureg.inch_H2O_39F # Based on the DP gage documentation
stopwatch = 0.01 * ureg.second # Photograph the stopwatch to achieve this
CO = 0.5 * 1.e-6 # Based on PEMS documentation
CO2 = 50 * 1.e-6 # Basedo on PEMS documentation
graduatedcyl = 0.5 * ureg.milliliter # Based on the resolution of the cylinder
ethtol = 0.03 * ureg.megajoule / ureg.kilogram # Based on article

"""
MEASUREMENT VALUES START HERE -------------------------------------------------
"""

# Add in all of the measured masses
fuelmass = 41 * ureg.gram # Measured value
potmeasurement = 99 * ureg.gram # Measured value
watermassmeasurement = 590 * ureg.gram # Measured value
ethanolmass = 274 * ureg.gram # Measured value
bottletotal = 304 * ureg.gram # Measured value
hotwatermass = 435 * ureg.gram # Measured value

# Poured volume of fuel in liters
fuelvolume = 50 * ureg.milliliter # Measured value

# Temperature measurements
starttemp = 298.45 * ureg.kelvin # Measured value
tempreading = Q_(75.4, ureg.degF) # Measured value
ambienttemp = tempreading.to(ureg.kelvin) # Convert to kelvin
finaltemp = 373.15 * ureg.kelvin # Solution was boiling

# Add in the pressure measurements of the sampling tube
# Will be input and output in inches of water
pressurediff = 0.456276 * ureg.inch_H2O_39F # Measured value

# Absolute pressure in Pascals
barpressure = 1024. * ureg.millibar # Measured value

# Exhaust gas constituent measurements
# Points 120 to 607 were measured
COlevel = 15.8902 * 1.e-6 * ureg.mol/ureg.mol # Measured value
CO2level = 283.2812 * 1.e-6 * ureg.mol/ureg.mol # Measured value

# Add in burn time, in seconds
minutes = 24. * ureg.minute # Measured value
seconds = 14. * ureg.second # Measured value

"""
CALCULATION STARTS HERE -------------------------------------------------------
"""

# Add uncertainties to all of the measured values
fuelv = fuelvolume.plus_minus(graduatedcyl)
methanol = ethanolmass.plus_minus(scale1)
mbottle = bottletotal.plus_minus(scale1)
mpot0 = potmeasurement.plus_minus(scale1)
mpot1 = watermassmeasurement.plus_minus(scale1)
mpot2 = hotwatermass.plus_minus(scale1)
deltaP = pressurediff.plus_minus(press)
Pabsolute = barpressure.plus_minus(barometer)
COppm = COlevel.plus_minus(CO)
CO2ppm = CO2level.plus_minus(CO2)
burntime = (minutes + seconds).plus_minus(1.).to(ureg.second)
T0 = starttemp.plus_minus(thermocouples)
Tambient = ambienttemp.plus_minus(thermometer)
TE = finaltemp.plus_minus(thermocouples)
Uethanol = Uethanol.plus_minus(ethtol)

# Plug in each of the equations and solve.  The bottom of the tree is computed
# first in order to allow computation of the higher level tolerances.
mfuel = (fuelv * ethrho90).to(ureg.kilogram)
rho_air = (Pabsolute / (R_air * Tambient)).to(ureg.kilogram / ureg.meter**3)
mwater1 = (mpot1 - mpot0).to(ureg.kilogram)
mwater2 = (mpot2 - mpot0).to(ureg.kilogram)
mevap = (mpot1 - mpot2).to(ureg.kilogram)
Cx = methanol / mbottle
Ux = (Uethanol * Cx).to(ureg.joule / ureg.kilogram)
Ereleased = (mfuel * Ux).to(ureg.joule)
EH2Oevap = (mevap * hv).to(ureg.joule)
deltaT = (TE - T0).to(ureg.kelvin)
N = (215. * (ureg.ft**3) /
     (ureg.minute*(ureg.inch_H2O_39F**0.5))).plus_minus(0)
Vdottube = (N*(deltaP**0.5)).to(ureg.meter**3 / ureg.second)
C = (1.5998468653396583 * ureg.meter/ureg.meter).plus_minus(
     0.02159726001744998) # Calibrated value
Vdot = Vdottube * C
V = (Vdot * burntime).to(ureg.meter**3)
COpartial = Pabsolute * COppm
CO2partial = Pabsolute * CO2ppm
nco = ((NA * COpartial * V) / (R * Tambient)).to('')
nco2 = ((NA * CO2partial * V) / (R * Tambient)).to('')
nethanol = NA * mfuel * Cx / Methanol
EH2Oheat = mwater2*cpw*deltaT + mpot0*cpm*deltaT
hc = (EH2Oheat + EH2Oevap) / Ereleased
epsilon = nco / (2 * nethanol)
epsilon2 = nco2 / (2 * nethanol)
sigma = ((nco * MCO / NA) / (EH2Oevap +
                             EH2Oheat)).to(ureg.gram / ureg.megajoule)

print('CONCENTRATION: {:.1f}%'.format(Cx.n * 100))
print('Sunday Value')
print('-----------------------------')
print('THERMAL EFFICIENCY (h_c)')
print('Nominal value: {:.2f}%'.format(hc.n * 100.))
print('Relative uncertainty: +/- {:.2f}%'.format(hc.rel *100.))
print('-----------------------------')
print('CO PRODUCTION')
print('Nominal value: {:.2f}%'.format(epsilon.n * 100.))
print('Relative uncertainty: +/- {:.2f}%'.format(epsilon.rel *100.))
print('-----------------------------')
print('CO2 PRODUCTION')
print('Nominal value: {:.2f}%'.format(epsilon2.n * 100.))
print('Relative uncertainty: +/- {:.2f}%'.format(epsilon2.rel *100.))
print('-----------------------------')
print('SPECIFIC EMMISSIONS')
print('Nominal value: {:.2f}'.format(sigma.n))
print('Relative uncertainty: +/- {:.2f}%'.format(sigma.rel *100.))

