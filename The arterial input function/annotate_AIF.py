import numpy as np
import os
from imagedata.series import Series
import matplotlib.pyplot as plt
from datetime import date

'''
This script is used for annotating an AIF using the imagedata Python library. 
'''

patient_path = '' #path to patientdata
dce_path = '' #path to DCE-MRI data
aif_path = '' #output path 

#Read DCE-MRI data 
DCE = Series(dce_path, 'time') 
timeline = DCE.timeline

#Subtract S0 
DCE = DCE.astype('float32')
S0 = np.mean(DCE[0:5,:,:,:], axis=0)
dynim = DCE.copy()
for k in range(len(timeline)):
    dynim[k,:,:,:] = DCE[k,:,:,:] - S0

#Chose AIF
AIF_roi = dynim.get_roi(follow=True)

#Plot AIF over time (NB: Always check that the AIF curve looks right)
AIF_np = np.array(AIF_roi, dtype='float')
non_zero_in_AIF = np.nonzero(AIF_np) #index of nonzero elements in mask image
AIF_vals = dynim[:, non_zero_in_AIF[1], non_zero_in_AIF[2], non_zero_in_AIF[3]]
AIF_vals = np.asarray(np.mean(AIF_vals, axis=1), dtype='float')
np.savetxt("H:/data/master endometrial data/198/AIF_198.csv", AIF_vals, delimiter=",")
plt.plot(timeline, AIF_vals)
plt.grid()
plt.title('Arterial input function')
plt.savefig(os.path.join(patient_path, 'AIF_plot.jpg'))
plt.show()

#Save AIF ROI as dicom
AIF_roi.header.seriesDescription = DCE.header.seriesDescription + ';AIF'
today = date.today()
d1 = today.strftime("%Y%m%d")
AIF_roi.tags = dynim.tags
AIF_roi.setDicomAttribute('AcquisitionDate', d1)
AIF_roi.setDicomAttribute('ContentDate', d1)
AIF_roi.setDicomAttribute('SeriesDate', d1)
AIF_roi.setDicomAttribute('StudyDate', d1)
AIF_roi.write(aif_path)

