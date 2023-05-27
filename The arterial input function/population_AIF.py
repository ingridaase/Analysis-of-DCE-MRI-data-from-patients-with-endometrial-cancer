import numpy as np
from numpy import savetxt
import os
from imagedata.series import Series
import matplotlib.pyplot as plt 

data_path = 'H:/data/master endometrial data/'
#patient_list = ['118', '122', '042', '136', '170', '255', '265', '198', '229', '135', '246']
patient_list = ['011', '036', '042', '086', '097', '118', '122', '135', '136', '150', '170', '191', '198','199', '203', '229', '246', '251', '255', '265']
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
    
#Mean of 11 AIF
mean_AIF = n_AIF.mean(axis=0)
print(mean_AIF)
savetxt('H:/data/master endometrial data/Pop_AIF/aif_data3.csv', mean_AIF, delimiter=',')

plt.plot(mean_AIF)
plt.title('Arterial Input Function (AIF)')
plt.xlabel('time [s]')
plt.ylabel('Signal intensity [a.u]')
plt.savefig('H:/data/master endometrial data/Pop_AIF/pop_aif_all.jpg')