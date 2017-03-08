#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 10:09:40 2017

@author: zaigerm

This python code is meant to read the data collected for the HotWire...
...measurments lab.

HOW TO USE>>>>>

1) Type in the filepath to the directory of what data you want to read
2) delete the .csv from the filenames list "del filenames[index of .csv]
    You can look in filenames prior to the del comment and make sure.

and if linked to a plotter you can plot the voltages for all of the lvm files

"""
#import lvm
import numpy as np
import scipy as sc
import os
from os.path import isfile, join
import matplotlib.pyplot as plt
#from graph import graph
#from lvm_reader import Lvm_read

def double_sort(l1, l2):
    a = l1
    a2 = l2
    #b = a[:]

    for j in range(len(l1)-1, 0, -1):
    #while True:
        #if is_sorted(a) == True:
            #break
        for i in range(len(a)-1):
            if a[i] > a[i+1]:
                holder1 = a[i]
                holder2 = a2[i]
                a[i] = a[i+1]
                a2[i] = a2[i+1]
                a[i+1] = holder1
                a2[i+1] = holder2

    return a, a2


def reader():
    filepath = '/nfs/stak/students/a/alfermaa/Desktop/Classes/ME552/Lab_2/Data/OHR_122/Calibration'
    filenames = [f for f in os.listdir(filepath) if isfile(join(filepath, f))]

    filenames = sorted(filenames)

    #del filenames[3] # This removes the .csv File from the filenames list since its
    # set up differently then the others.

    aa = np.arange(11)
    R0 = 5.85 # resistance at zero
    RcR0 = 1.33 # Difference between zero and 100C  resistance
    Rop = 7.44 # Operation Resistance
    Rint = 0.5 #given value
    Rw = Rop + Rint # Resistance of the wire


    Tw = 100 * (Rop - R0 / RcR0)
    Tatm = 20.3 +273.15

    T_a = np.array([20.5, 20.8, 20.7, 20.6, 20.6, 20.5, 20.6, 20.6, 20.6, 20.6, 20.6]) +273.15
    DP122 = (np.array([0.0032, 0.0086 ,0.0648, 0.0038, 0.0528, 0.0054, 0.0045, 0.0184, 0.007, 0.0034, 0.0334]) - 0.0032)*1000
    DPatm = 0.0032*1000

    DP122 = np.array(sorted(DP122))
    #T = sorted(DP122)

    R = 286.9# sc.constants.gas_constant
    utm = (2*DPatm*R*Tatm/101325.0)**0.5
    ucal =(2*DP122*R*T_a/101325.0)**0.5

    avgV = []
    stdv = []
    keys = np.arange(11)

    for name in filenames:
        voltage = []
        with open(filepath +'/' + name, 'r') as f:
            line_15 = f.readlines()[23:]
            for line in line_15:
                lines = line.strip()
                voltage.append(float(lines))
        vavg = np.average(voltage)
        avgV.append(vavg)
        stds = np.std(voltage)
        stdv.append(stds)
        f.close()

    A = avgV[0] / (Rw *( Tw -Tatm))
    B = (avgV / (Rw * (Tw -T_a) - A))/(ucal**3)
    #print(B)
    #print ucal

    u = (((avgV / (Rw * (Tw -T_a) - A))/(B))**0.3333333)
    # print(avgV)



#    a = double_sort(np.absolute(avgV), stdv)
#    ps = double_sort(np.absolute(avgV), DP122)
#    ts = double_sort(np.absolute(avgV), T_a)
    DP122 = sorted(DP122)

#    volts = np.array(a[0])
#    std2 = np.array(a[1]) * 2.0
#    ts = np.array(ts[1])
#    ps = np.array(ps[1])
    #print DP122
    #print ps
    #print a
    #print sorteds
    #avgV = sorted(np.absolute(avgV))
    u = sorted(u)
    ucal =sorted(ucal)

    volts = abs(np.array(avgV))

    ustds = np.array(stdv) * 2.

    return volts, ustds, T_a, DP122

if __name__ == '__main__':
    V, u_stds, Tambreadings, Pressreadings = reader()

#    plt.figure()
#    plt.errorbar(a[0], u, yerr = std2, fmt ='o')
#    z = np.polyfit(a[0], u, 3)
#    p = np.poly1d(z)
#    plt.plot(a[0], p(a[0]), 'r--')
#    #plt.plot(a[0], 'ro' )
#    plt.show()    #graph()
#    print "y=%.6fx^3+(%.6f)x^2+(%.6f)x+(%.6f)"%(z[0],z[1], z[2], z[3])
