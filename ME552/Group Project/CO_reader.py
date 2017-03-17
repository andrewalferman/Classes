#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 13:38:45 2017

@author: andrewalferman
"""

"""
Reader file for the ME552 Winter 2017 Group 1 Term Project
This file grabs all of the raw data from the csv files that the portable
emissions measurement system (PEMS) software generates, and processes it into a
more user friendly file.  It also generates plots of all of the data, if
desired.

The file is currently configured to not process some of the data from the files
due to missing data from the PEMS output.

"""

import os as os
import numpy as np
import csv as csv
import pylab as pyl
import pint

def concentrations(CO,CO2,**mod):
    if ('startpoint' in mod):
        CO = CO[mod['startpoint']:]
        CO2 = CO2[mod['startpoint']:]
    COtotal = 0.
    CO2total = 0.
    for i in range(len(CO)):
        COtotal += float(CO[i])
        CO2total += float(CO2[i])
    COtotal = (COtotal / len(CO)) - float(CO[0])
    CO2total = (CO2total / len(CO2)) - float(CO2[0])
    return COtotal, CO2total

def avgval(vals):
    total = 0.
    for i in range(len(vals)):
        total += float(vals[i])
    return total/len(vals)

"""
-------------------------------------------------------------------------------
If all goes well, you should only need to update the portion of the code
between these two commented lines.
"""
# Enter the filename that you want to parse. Don't add in the ".csv"
filename = '90_11mar2017_1401'

# Make this a 1 if you want to save the data or plot the data
savefile = 0
plotfigs = 0

"""
-------------------------------------------------------------------------------
"""
# Get the current path to add to the filename
currentpath = os.getcwd()

# Need to look at each file and determine which row has the multipliers, and
# which row is the last row of column headers.
multiplierrows = {'70_12mar2017_1': 4, '80_12mar2017_1': 5,
                  '90_12mar2017_1': 5, '100_12mar2017_1': 4,
                  '70_11mar2017_1530': 4, '80_11mar2017_1440': 5,
                  '90_11mar2017_1401': 4
                  }

# Use this dictionary for all of the starting point values for emmission avg
starts = {'70_12mar2017_1': 0, '80_12mar2017_1': 40, '90_12mar2017_1': 120,
          '100_12mar2017_1': 140, '70_11mar2017_1530': 80,
          '80_11mar2017_1440': 159, '90_11mar2017_1401': 40
          }

# Initialize lists to append the corrected data to
sec = []
CO = []
CO2 = []
PMred = []
flow = []
fluetemp = []
ThermoCo = []
PMgreen = []
FlowGrav = []
COtemp = []
gastemp = []
rlas = []
temp = []
rboxtemp = []
glastemp = []
gboxtemp = []
#gasboxRH = []
#PMgrenRH = []
chiptemp = []

# Read the data
with open(currentpath + '/Data/' + filename + '.csv', 'r', newline='') as f:
    reader = csv.reader(f, delimiter=',')
    # Iterate accross the rows, and use a counter to grab the information
    rownum = 1
    for row in reader:
        # Grab the multipliers
        if rownum == multiplierrows[filename]:
            COmult = float(row[1])
            CO2mult = float(row[2])
            PMredmult = float(row[3])
            flowmult = float(row[4])
            fluetempmult = float(row[5])
            ThermoComult = float(row[6])
            PMgreenmult = float(row[7])
            FlowGravmult = float(row[8])
            COtempmult = float(row[9])
            gastempmult = float(row[10])
            rlasmult = float(row[11])
            tempmult = float(row[12])
            rboxtempmult = float(row[13])
            glastempmult = float(row[14])
            gboxtempmult = float(row[15])
            #gasboxRHmult = float(row[16])
        # This if statement cuts off all the info at the top of the file
        if rownum > multiplierrows[filename] + 3:
            # Convert all of the data from strings to floats
            row = [float(e) for e in row if e]
            # Test to see if the row contains most of the useful information
            try:
                test = row[8]
                realrow = True
            except IndexError:
                realrow = False
            # Append all of the useful data from the larger rows
            if realrow == True:
                sec.append(row[0])
                CO.append(row[1])
                CO2.append(row[2])
                PMred.append(row[3])
                flow.append(row[4])
                fluetemp.append(row[5])
                ThermoCo.append(row[6])
                PMgreen.append(row[7])
                FlowGrav.append(row[8])
                COtemp.append(row[9])
                gastemp.append(row[10])
                rlas.append(row[11])
                temp.append(row[12])
                rboxtemp.append(row[13])
                glastemp.append(row[14])
                gboxtemp.append(row[15])
                #gasboxRH.append(row[16])
            # Append the data from the smaller rows
            else:
                #PMgrenRH.append(row[-2])
                chiptemp.append(row[-1])
        else:
            pass
        # Be sure to progress the counter
        rownum += 1

# Incorporate the multiplier on each of the lists, as applicable
CO = [x * COmult for x in CO]
CO2 = [x * CO2mult for x in CO2]
PMred = [x * PMredmult for x in PMred]
flow = [(x - 10550.)/41000. for x in flow]
fluetemp = [x * fluetempmult for x in fluetemp]
ThermoCo = [x * ThermoComult for x in ThermoCo]
PMgreen = [x * PMgreenmult for x in PMgreen]
FlowGrav = [x * FlowGravmult for x in FlowGrav]
COtemp = [x * COtempmult for x in COtemp]
gastemp = [x * gastempmult for x in gastemp]
rlas = [x * rlasmult for x in rlas]
temp = [x * tempmult for x in temp]
rboxtemp = [x * rboxtempmult for x in rboxtemp]
glastemp = [x * glastempmult for x in glastemp]
gboxtemp = [x * gboxtempmult for x in gboxtemp]
#gasboxRH = [x * gasboxRHmult for x in gasboxRH]

# Put headers on each of the rows (which will be converted to columns)
sec = np.hstack(('seconds',sec))
CO = np.hstack(('CO',CO))
CO2 = np.hstack(('CO2',CO2))
PMred = np.hstack(('PMred',PMred))
flow = np.hstack(('DP',flow))
fluetemp = np.hstack(('fluetemp',fluetemp))
ThermoCo = np.hstack(('ThermoCo',ThermoCo))
PMgreen = np.hstack(('PMgreen',PMgreen))
FlowGrav = np.hstack(('FlowGrav',FlowGrav))
COtemp = np.hstack(('COtemp',COtemp))
gastemp = np.hstack(('gastemp',gastemp))
rlas = np.hstack(('rlas',rlas))
temp = np.hstack(('temp',temp))
rboxtemp = np.hstack(('rboxtemp',rboxtemp))
glastemp = np.hstack(('glastemp',glastemp))
gboxtemp = np.hstack(('gboxtemp',gboxtemp))
#gasboxRH = np.hstack(('gasboxRH',gasboxRH))
#PMgrenRH = np.hstack(('PMgrenRH',PMgrenRH))
chiptemp = np.hstack(('chiptemp',chiptemp))

# Putting it like this so that we can easily forego columns
A = np.array([sec,
              CO,
              CO2,
              PMred,
              flow,
              fluetemp,
              ThermoCo,
              PMgreen,
              FlowGrav,
              COtemp,
              gastemp,
              rlas,
              temp,
              rboxtemp,
              glastemp,
              gboxtemp,
              #gasboxRH,
              #PMgrenRH,
              chiptemp
              ])
# Need to transpose the data for it to be written in a csv, if desired
A = np.transpose(A)

# Find the region of interest based on where we have values for ThermoCo
imin = 9999
imax = 0
for i in range(1,len(ThermoCo)):
    # 6553.6 chosen because that's the value that is returned when the
    # thermocouple is disconnected
    if float(ThermoCo[i]) != 6553.6 and i < imin:
        imin = i
    if float(ThermoCo[i]) != 6553.6 and i > imax:
        imax = i

# Save the data in a csv file, if desired
if savefile == 1:
    with open(filename + '_corrected_flow.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(A)

# Plot the data, if desired
if plotfigs == 1:
    # Plot only the data we want
    xmin = float(sec[imin])
    xmax = float(sec[imax])
    size = (6, 4.5)
    fs = 14
    dots = 900
    # Start at 1 because seconds vs. seconds is useless
    for i in range(1,len(A[0,:])):
        pyl.figure(i-1, figsize=size, dpi=dots)
        pyl.plot(sec[imin:imax], A[imin:imax,i])
        pyl.xlabel('Time (Seconds)', fontsize=fs)
        pyl.grid(True)
        pyl.xlim(xmin,xmax)

    # Add the figure specific labels and settings
    # Keep in mind that it is not configured to plot gasboxRH or PMgrenRH
    # currently
    pyl.figure(0)
    pyl.ylabel('CO Concentration (ppm)', fontsize=fs)
    pyl.ylim(0,)
    pyl.figure(1)
    pyl.ylabel('CO2 Concentration (ppm)', fontsize=fs)
    pyl.figure(2)
    pyl.ylabel('PMred Value', fontsize=fs)
    pyl.figure(3)
    pyl.ylabel('DP Value', fontsize=fs)
    pyl.figure(4)
    pyl.ylabel('Flue Temp Value', fontsize=fs)
    pyl.figure(5)
    pyl.ylabel('Pot Temperature ($^\circ$C)', fontsize=fs)
    pyl.figure(6)
    pyl.ylabel('PMgreen Value', fontsize=fs)
    pyl.figure(7)
    pyl.ylabel('FlowGrav Value', fontsize=fs)
    pyl.figure(8)
    pyl.ylabel('CO Temp Value', fontsize=fs)
    pyl.figure(9)
    pyl.ylabel('Gas Temp Value', fontsize=fs)
    pyl.figure(10)
    pyl.ylabel('RLAS Value', fontsize=fs)
    pyl.figure(11)
    pyl.ylabel('Temp Value', fontsize=fs)
    pyl.figure(12)
    pyl.ylabel('Rbox Temp Value', fontsize=fs)
    pyl.figure(13)
    pyl.ylabel('GLAS Temp Value', fontsize=fs)
    pyl.figure(14)
    pyl.ylabel('Gbox Temp Value', fontsize=fs)
    pyl.figure(15)
    pyl.ylabel('Chip Temp Value', fontsize=fs)

    # Show the goods
    pyl.show()

# Display the average CO and CO2 concentrations as well as the average DP
COavg, CO2avg = concentrations(CO[imin:imax],
                               CO2[imin:imax],
                               startpoint=(starts[filename]-imin + 1)
                               )

# Some fluctuations in one file that we want to fix
if filename == '100_12mar2017_1':
    imin += 4

# Find the temperature of the water at the start and stop of the WBT
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
Tempmin = Q_(float(ThermoCo[imin]), ureg.degC)
Tempmax = Q_(float(ThermoCo[imax]), ureg.degC)

# Find the time that it takes to boil
boilthreshold = 98.0
for i in range(starts[filename]+1, len(ThermoCo[imin:imax])):
    if float(ThermoCo[i]) >= boilthreshold:
        boilnum = i
        break
    elif i == imax - 1:
        boilnum = False
if boilnum != False:
    boiltime = float(sec[boilnum]) - float(sec[starts[filename]+1])

print(filename)
print('First point: {}'.format(starts[filename]))
print('Last point: {}'.format(imax))
print('CO Average: {:.4f} ppm'.format(COavg))
print('CO2 Average: {:.4f} ppm'.format(CO2avg))
print('DP Average: {:.6f} inH2O'.format(avgval(flow[1:])))
print('Start Temp: ' + str(Tempmin.to(ureg.kelvin)))
print('Final Temp: ' + str(Tempmax.to(ureg.kelvin)))
if boilnum != False:
    print('Boil time: {} seconds'.format(boiltime))
else:
    print('Water did not reach {} degC'.format(boilthreshold))
