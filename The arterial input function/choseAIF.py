import numpy as np
import os
from imagedata.series import Series
import matplotlib.pyplot as plt
from datetime import date

#Select a AIF ROI from the dynamic image 
'''
#Pas 011
patient_path = 'H:/data/master endometrial data/011/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 036
patient_path = 'H:/data/master endometrial data/036/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')
#Pas 042
patient_path = 'H:/data/master endometrial data/042/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 086
patient_path = 'H:/data/master endometrial data/086/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 097
patient_path = 'H:/data/master endometrial data/097/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Patient 118
patient_path = 'H:/data/master endometrial data/118/'
aif_path = os.path.join(patient_path, 'AIF')
dce_path = os.path.join(patient_path, 'dce')
mask_path = os.path.join(patient_path, 'mask')

#Patient 122
patient_path = 'H:/data/master endometrial data/122/'
dce_path = os.path.join(patient_path, 'dce')
vibe_path = os.path.join(patient_path, 'vibe')
mask_path = os.path.join(patient_path, 'mask')
aif_path = os.path.join(patient_path, 'AIF')

#Patient 135
patient_path = 'H:/data/master endometrial data/135/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
vibe_path = os.path.join(patient_path, 'vibe')
aif_path = os.path.join(patient_path, 'AIF')

#Patient 136
patient_path = 'H:/data/master endometrial data/136/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
vibe_path = os.path.join(patient_path, 'vibe')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 150
patient_path = 'H:/data/master endometrial data/150/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Patient 170
patient_path = 'H:/data/master endometrial data/170/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Patient 191
patient_path = 'H:/data/master endometrial data/191/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 198
patient_path = 'H:/data/master endometrial data/198/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 199
patient_path = 'H:/data/master endometrial data/199/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 203
patient_path = 'H:/data/master endometrial data/203/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 229
patient_path = 'H:/data/master endometrial data/229/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 246
patient_path = 'H:/data/master endometrial data/246/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 251
patient_path = 'H:/data/master endometrial data/251/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')

#Pas 255
patient_path = 'H:/data/master endometrial data/255/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')
'''
#Pas 265
patient_path = 'H:/data/master endometrial data/265/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')
'''
#Pas 387
patient_path = 'H:/data/master endometrial data/387/'
mask_path = os.path.join(patient_path, 'mask')
dce_path = os.path.join(patient_path, 'dce')
aif_path = os.path.join(patient_path, 'AIF')
'''
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

'''
#If plotting, then store the time and slice for the image the AIF was chosen from 
print('Enter time for AIF:')
AIF_time = input()
AIF_time = int(AIF_time)
print('Enter slice for AIF:')
AIF_slice = input()
AIF_slice = int(AIF_slice)

AIF_np = np.array(AIF_roi, dtype='float')
#Plot chosen AIF on dce image
fig = plt.figure()
plt.imshow(dynim[AIF_time, AIF_slice,:,:], cmap='gray')
plt.imshow(AIF_np[AIF_time, AIF_slice,:,:], alpha=0.3)
plt.title('AIF ROI')
plt.axis('off')
plt.show()
'''
