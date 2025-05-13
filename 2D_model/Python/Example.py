
import numpy as np
import scipy.io as sc
import matplotlib.pylab as plt
from URB import URB

#%% Load of experimental results and model parameters

Pmt=np.loadtxt('../Data/Parameters.txt')
variable=['a1','a2','a3','fs1','ps1','fs2','ps2','fs3','ps3','fm','pm','uy','fy','B','n','fmax','pphi']

var=dict()
i = 0
for vnam in variable:
    var[vnam]=float(Pmt[i])
    i += 1

data = sc.loadmat('../Data/Test.mat')
var['h']= float(data['h'][0][0])
d=data['d'][0]
f=data['f'][0]


#%% Solution

f_model=list()
for i in d:
    (ftemp,var)=URB(i,var)
    if (ftemp == 'Failed'):
        break
    f_model.append(ftemp)


#%% Figure
if (ftemp != 'Failed'):
    plt.figure
    plt.plot(d,f,color='k',label=r'Experimental')
    plt.plot(d,f_model,color='r',label=r'Proposed model')
    plt.ylim(-45.0, 45.0)
    plt.xlim(-25.0, 25.0)
    plt.legend()
    plt.grid(True)
    plt.title('Comparison between numerical and experimental results')
    plt.xlabel('Deformation [cm]')
    plt.ylabel('Force [tonf]')

plt.savefig('../Figures/Python.pdf')