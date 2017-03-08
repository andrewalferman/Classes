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
import os
from os.path import isfile, join
import matplotlib.pyplot as plt
#from lvm_reader import Lvm_read

filepath = '/nfs/stak/students/a/alfermaa/Desktop/Classes/ME552/Lab_2/Data/OHR_122/Calibration'
filenames = [f for f in os.listdir(filepath) if isfile(join(filepath, f))]

del filenames[3]




for name in filenames:
    voltage = []
    with open(filepath +'/' + filenames[0], 'r') as f:
        line_15 = f.readlines()[23:]
        for line in line_15:
            lines = line.strip()
            voltage.append(float(lines))

print(min(voltage))