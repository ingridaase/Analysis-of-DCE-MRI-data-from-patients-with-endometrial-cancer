import numpy as np
from numpy import savetxt
import os
from imagedata.series import Series
import matplotlib.pyplot as plt 

data_path = '' #Path to AIF data

patient_list = [] #List of patients to retrieve AIF data from 
n_patients = len(patient_list)
n_AIF = np.zeros((n_patients, 160))

AIF_series = []
dce_series = []
for index,item in enumerate(patient_list): 
    patient_path = os.path.join(data_path, item)
    aif_path = os.path.join(patient_path, 'AIF')
    dce_path = os.path.join(patient_path, 'dce')
    
    AIF = Series(aif_path, 'time')
    AIF_series.append(AIF)
    AIF_np = np.asarray(AIF)
    DCE = Series(dce_path, 'time')
    timeline = DCE.timeline
    
    #Subtract pre-contrast signal from signal
    DCE = DCE.astype('float32')
    #Subtract pre-contrast signal from signal
    S0 = np.mean(DCE[0:5,:,:,:], axis=0)
    dynim = DCE.copy()
    for k in range(len(timeline)):
        dynim[k,:,:,:] = DCE[k,:,:,:] - S0
    
    #Find mean values inside AIF mask 
    non_zero_in_AIF = np.nonzero(AIF_np) #index of nonzero elements in mask image
    AIF_vals = dynim[:, non_zero_in_AIF[1], non_zero_in_AIF[2], non_zero_in_AIF[3]]
    AIF_vals = np.asarray(np.mean(AIF_vals, axis=1), dtype='float')
    n_AIF[index,:] = AIF_vals
    
#Mean of 20 AIF
mean_AIF = n_AIF.mean(axis=0)
print(mean_AIF)
savetxt('', mean_AIF, delimiter=',') #Choose path where AIF should be saved

plt.plot(mean_AIF)
plt.title('Arterial Input Function (AIF)')
plt.xlabel('time [s]')
plt.ylabel('Signal intensity [a.u]')
plt.savefig('') #Choose path where figure should be saved