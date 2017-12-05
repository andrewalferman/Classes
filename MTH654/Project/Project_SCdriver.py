#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 15:02:15 2016

@author: mattz

Heavily modified by Andrew Alferman for MTH 654 term project
12/1/17
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import os as os
import csv as csv


def loadpasrdata(num):
    """Load the initial conditions from the PaSR files."""
    pasrarrays = []
    for i in range(num):
        filepath = os.path.join(os.getcwd(),
                                'pasr_out_h2-co_' +
                                str(i) +
                                '.npy')
        filearray = np.load(filepath)
        pasrarrays.append(filearray)
    return np.concatenate(pasrarrays, 1)


def getH2O2ics(pasr, p, t):
    """Unpack the initial conditions in a legible format."""
    temperature = pasr[t, p, :][1]
    pressure = pasr[t, p, :][2]
    H = pasr[t, p, :][3]
    H2 = pasr[t, p, :][4]
    O = pasr[t, p, :][5]
    OH = pasr[t, p, :][6]
    H2O = pasr[t, p, :][7]
    O2 = pasr[t, p, :][8]
    HO2 = pasr[t, p, :][9]
    H2O2 = pasr[t, p, :][10]
    N2 = pasr[t, p, :][11]
    AR = pasr[t, p, :][12]
    HE = pasr[t, p, :][13]
    CO = pasr[t, p, :][14]
    CO2 = pasr[t, p, :][15]
    stddevvector = [0.1 * i for i in pasr[t, p, 3:11]]
    return temperature, pressure, H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, \
        HE, CO, CO2, stddevvector


def reactionsetup(time, pasr, p, t, randomization):
    """Set up the reaction."""
    # reset the ignition time flag and temperature
    T_burn = np.zeros(len(time))
    # utilize a constant pressure reactor from cantera
    # to set equilivalence ratios.
    gas = ct.Solution('chem.cti')
    # gas.TPX = 1000.0, ct.one_atm, ...
    #   'CH4:{0},O2:2.0,N2:7.52'.format(phi[eq])
    temperature, pressure, H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, HE, \
        CO, CO2, stddevvector = getH2O2ics(pasr, p, t)
    if randomization:
        H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, HE, CO, CO2, randomvector = \
            randomizedH2O2massfractions(H, H2, O, OH, H2O, O2, HO2,
                                        H2O2, N2, AR, HE, CO, CO2)
    canterastring = ('H:{},H2:{},O:{},OH:{},H2O:{},O2:{},' +
                     'HO2:{},H2O2:{},N2:{},AR:{},HE:{},CO:{},' +
                     'CO2:{}').format(
                     H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, HE,
                     CO, CO2)
    gas.TPX = temperature, ct.one_atm, canterastring
    react = ct.IdealGasConstPressureReactor(gas)
    soln = ct.ReactorNet([react])
    return T_burn, react, soln, stddevvector


def randomizedH2O2massfractions(H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, HE,
                                CO, CO2):
    """Randomize the state."""
    physical = False
    while not physical:
        randomvector = np.random.normal(size=8)
        stddevfraction = 0.1
        H = H + H * stddevfraction * randomvector[0]
        H2 = H2 + H2 * stddevfraction * randomvector[1]
        O = O + O * stddevfraction * randomvector[2]
        OH = OH + OH * stddevfraction * randomvector[3]
        H2O = H2O + H2O * stddevfraction * randomvector[4]
        O2 = O2 + O2 * stddevfraction * randomvector[5]
        HO2 = HO2 + HO2 * stddevfraction * randomvector[6]
        H2O2 = H2O2 + H2O2 * stddevfraction * randomvector[7]
        if ((H > 0.0 and H < 1.0) and
           (H2 > 0.0 and H2 < 1.0) and
           (O > 0.0 and O < 1.0) and
           (OH > 0.0 and OH < 1.0) and
           (H2O > 0.0 and H2O < 1.0) and
           (O2 > 0.0 and O2 < 1.0) and
           (HO2 > 0.0 and HO2 < 1.0) and
           (H2O2 > 0.0 and H2O2 < 1.0)):
                physical = True
    N2 = 1 - (H + H2 + O + OH + H2O + O2 + HO2 + H2O2)
    return H, H2, O, OH, H2O, O2, HO2, H2O2, N2, AR, HE, \
        CO, CO2, randomvector


def gaussianpdf(x, mu, sigma):
    """Return the probability density function.

    Inputs:
    x - sample value
    mu - average value
    sigma - standard deviation
    """
    return np.exp(-0.5*((x - mu)/sigma)**2)/(sigma * np.sqrt(2*np.pi))


t = 1.0
dt = 0.001
n = int(t/dt)
# initalize time and a storage for when the temperature spikes
time = np.linspace(0, t, n)
change = 0.0
# plt.figure(num=None, figsize=(12, 8), dpi=900, facecolor='w', edgecolor='k')

pasr = loadpasrdata(1)
numparticles = len(pasr[0, :2, 0])
numtsteps = len(pasr[:, 0, 0])

# Forget what we just did and select the best particle for demonstration
# numparticles = [83]
# numtsteps = [881]

# ignited = False

"""Get all of the sparse herm weights and points
# Data obtained by running the following Matlab code:

NSP = 8;
level = 2;
point_num = sparse_grid_herm_size(NSP, level);
[xw, xs] = sparse_grid_herm(NSP, level, point_num);

csvwrite('weights_sparse_herm.csv', xw)
csvwrite('points_sparse_herm.csv', xs)
"""
with open('weights_sparse_herm.csv') as csvfile:
    rows = csv.reader(csvfile, delimiter=',')
    for row in rows:
        weights = [float(i) for i in row]
points = []
with open('points_sparse_herm.csv') as csvfile:
    rows = csv.reader(csvfile, delimiter=',')
    for row in rows:
        points.append([float(i) for i in row])
points = np.array(points)

# Run the nominal reaction start to finish
T_burn, react, soln, stddevvector = reactionsetup(time, pasr, 83, 881, False)
for i in range(len(time)):
    soln.advance(time[i])
    # now for each time step pull store the temperature
    T_burn[i] = react.T
# Apply stochastic collocation
T_cases = []
for x in range(len(weights)):
    T_burnR, react, soln, stddevvectorR = \
        reactionsetup(time, pasr, 83, 881, True)
    for i in range(len(time)):
        soln.advance(time[i])
        # now for each time step pull store the temperature
        T_burnR[i] = react.T

plt.figure(0)
plt.plot(time, T_burn)
plt.xlabel('Time (sec)')
plt.ylabel('Temperature (K)')
plt.title('Nominal Case Time vs. Temperature')
plt.grid(b=True, which='both')
plt.savefig('nominalcase.png', dpi=600)

plt.show()
