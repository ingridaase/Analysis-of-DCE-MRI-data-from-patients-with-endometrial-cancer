from runtofts import *
from readdata import *
import numpy as np
import os
import json
from numpy import loadtxt


'''
Script used for running ETM locally for a population-based AIF. The changes needed to run ETM for manually annotated AIFs or for automatic AIFs 
are included in comments
'''

#Model options 
read_option = 'mask2dce'
model_option = 'voxelwise'

# load aif array i population-based AIF is used 
aif = loadtxt('H:/data/master/Pop_AIF/aif_data.csv', delimiter=',') #Load population AIF

data_path = '' #path to patient data
patient_list = ['1', '2', '3', '4', '5', '6'] #only example of a patientlist

#Run tofts method for all patients 
params = []
#AIFcost = [] If automatic AIF is used 
for i in patient_list: 
    patient_path = os.path.join(data_path, i)
    #aif_path = os.path.join(patient_path, 'AIF') If manually annotated AIFs are used
    mask_path = os.path.join(patient_path, 'mask')
    dce_path = os.path.join(patient_path, 'dce')
    opt_params = runtofts_popAIF(dce_path, mask_path, aif, read_option, model_option)
    #opt_params = runtofts_manualAIF(dce_path, mask_path, aif_path, read_option, model_option) If manually annotated AIFs are used
    
    #try: 
    #    opt_params, AIF_params = runtofts_autoAIF(dce_path, mask_path, read_option, model_option) if automatic AIF is used 
    #    params.append(opt_params)
    #    AIFcost.append(AIF_params)
    #except UnboundLocalError: 
    #    continue 
    params.append(opt_params)

# Writing to json
if model_option == 'average': 
    with open("H:/data/Results/PopAIF/output_average.json", "w") as outfile: 
        outfile.write(json.dumps(params, indent=4))
elif model_option == 'voxelwise':
    with open("H:/data/Results/PopAIF/output_voxelwise.json", "w") as outfile: 
        outfile.write(json.dumps(params, indent=4))