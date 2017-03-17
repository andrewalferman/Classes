#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:43:23 2017

@author: andrewalferman
"""

import pylab as pyl
import numpy as np

concen = [70.1, 80.3, 90.1, 100]
concenunc = [0.427, 0.403, 0.491, 0]

thermeff = [43.72, 43.55, 42.21, 38.96]
thermunc = [1.87, 1.70, 1.67, 1.56]

COproduc = [6.87, 6.84, 6.83, 5.02]
COuncert = [5.98, 5.55, 5.44, 3.67]

COperMJ = [6.45, 6.45, 6.64, 5.29]
COMJunc = [6.07, 5.60, 5.47, 3.72]

boiltime = [1052, 839, 686, 504]
boiltol = [3, 3, 3, 3]

fs = 14

pyl.figure(0, figsize=(8,6), dpi=600)
pyl.xlabel('Ethanol Concentration (%)', fontsize=fs)
pyl.ylabel('Thermal Efficiency (%)', fontsize=fs)
pyl.grid(True)
pyl.xlim(65,105)
pyl.errorbar(concen, thermeff,
             yerr=[thermeff[i]*thermunc[i]/100 for i in range(len(thermeff))],
             xerr=[concen[i]*concenunc[i]/100 for i in range(len(concen))],
             fmt='o'
             )
pyl.savefig('thermalefficiency.png')

pyl.figure(1, figsize=(8,6), dpi=600)
pyl.xlabel('Ethanol Concentration (%)', fontsize=fs)
pyl.ylabel('CO Production (%)', fontsize=fs)
pyl.grid(True)
pyl.xlim(65,105)
pyl.errorbar(concen, COproduc,
             yerr=[COproduc[i]*COuncert[i]/100 for i in range(len(COproduc))],
             xerr=[concen[i]*concenunc[i]/100 for i in range(len(concen))],
             fmt='o'
             )
pyl.savefig('coproduction.png')

pyl.figure(2, figsize=(8,6), dpi=600)
pyl.xlabel('Ethanol Concentration (%)', fontsize=fs)
pyl.ylabel('Specific CO Emmissions (g/MJ)', fontsize=fs)
pyl.grid(True)
pyl.xlim(65,105)
pyl.errorbar(concen, COperMJ,
             yerr=[COperMJ[i]*COMJunc[i]/100 for i in range(len(COperMJ))],
             xerr=[concen[i]*concenunc[i]/100 for i in range(len(concen))],
             fmt='o'
             )
pyl.savefig('specificemmissions.png')

pyl.figure(3, figsize=(8,6), dpi=600)
pyl.xlabel('Ethanol Concentration (%)', fontsize=fs)
pyl.ylabel('Boil Time (s)', fontsize=fs)
pyl.grid(True)
pyl.xlim(65,105)
pyl.errorbar(concen, boiltime,
             yerr=boiltol,
             xerr=[concen[i]*concenunc[i]/100 for i in range(len(concen))],
             fmt='o'
             )
pyl.savefig('boiltime.png')