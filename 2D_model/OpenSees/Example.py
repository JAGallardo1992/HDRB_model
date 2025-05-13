# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:14:08 2022

@author: jogallardo
"""

import subprocess
import pandas as pd
import matplotlib.pyplot as plt

#%% Function to load data

def read_response_2D(dir_files,case):
    
    dir_displ = dir_files+'/Displ_'+case+'.txt'
    dir_force = dir_files+'/Force_'+case+'.txt'
    
    displ = pd.read_csv(dir_displ, sep = ' ', names = ['Ind', 'd11', 'd12', 'd13', 'd21', 'd22', 'd23'])
    force = pd.read_csv(dir_force, sep = ' ', names = ['Ind', 'd11', 'd12', 'd13', 'd21', 'd22', 'd23'])
    displ.pop('Ind')
    force.pop('Ind')
    return displ, force

#%% Run the analysis

print('Runing HDRB...')
subprocess.run(['OpenSees','Run.tcl'], stderr=False)

#%% Plot the results

dpls_HDRB_D1, frc_HDRB_D1 = read_response_2D('Results','HDRB_D1')

plt.figure()
plt.plot(-dpls_HDRB_D1['d22'],-frc_HDRB_D1['d22'])
plt.grid(True)
plt.title('HDRB')
plt.show()












