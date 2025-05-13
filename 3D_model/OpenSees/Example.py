# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:14:08 2022

@author: jogallardo
"""

import subprocess
import pandas as pd
import matplotlib.pyplot as plt


#%% load data

def read_response_3D(dir_files,case):
    
    dir_displ = dir_files+'/Displ_'+case+'.txt'
    dir_force = dir_files+'/Force_'+case+'.txt'
    
    displ = pd.read_csv(dir_displ, sep = ' ', names = ['Ind', 'd11', 'd12', 'd13', 'd14', 'd15', 'd16','d21', 'd22', 'd23', 'd24', 'd25', 'd26',])
    force = pd.read_csv(dir_force, sep = ' ', names = ['Ind', 'd11', 'd12', 'd13', 'd14', 'd15', 'd16','d21', 'd22', 'd23', 'd24', 'd25', 'd26',])
    displ.pop('Ind')
    force.pop('Ind')
    return displ, force


#%% Run Analysis

print('Runing TSM...')
subprocess.run(['OpenSees','Run.tcl'], stderr=False)
dpls_TSM, frc_TSM = read_response_3D('Results','TSM')

#%% Plot results

plt.figure()
plt.plot(dpls_TSM['d22'],frc_TSM['d22'], label = 'X')
plt.plot(dpls_TSM['d23'],frc_TSM['d23'], '--', label = 'Y')
plt.grid(True)
plt.legend()
plt.xlabel('Shear deformation [cm]')
plt.ylabel('Shear force [tonf]')
plt.show()

plt.figure()
plt.plot(dpls_TSM['d21'],dpls_TSM['d22'])
plt.grid(True)
plt.xlabel('Axial deformation [cm]')
plt.ylabel('Axial force [tonf]')
plt.show()
