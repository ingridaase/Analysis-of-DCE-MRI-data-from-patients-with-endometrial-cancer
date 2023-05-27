#Import packages 
import numpy as np
import sys
import os
import json
from numpy import loadtxt
from imagedata import Study


#Import functions 
from runtofts import *

datafolder = sys.argv[1]

description = {}
with open(os.path.join(datafolder, "descr.json")) as f:
    description = json.load(f)[0]

output=os.path.join(datafolder,"output")
if not(os.path.exists(output)):
    try:
        os.mkdir(output,0o777)
    except OSError as error:
        (error)

patientID = description['PatientID']
    
aif = loadtxt('aif.csv', delimiter=',') #Population-based AIF

study = Study(datafolder+'/input', opts={'strict_values': False})
for uid in study: 
    series = study[uid]
    #if 'AIF' in series.seriesDescription: 
    #    aif = series
    if '12' in series.seriesDescription and 'AIF' not in series.seriesDescription:
        dce = series
    if 'JAD' in series.seriesDescription:
        mask_jad = series
    if 'KWL' in series.seriesDescription: 
        mask_kwl = series 



#Read option can be 'mask2dce' or 'dce2mask'
read_option = 'mask2dce'

model_option = 'average'

redcap_result = []

try: 
    opt_params = runtofts_popAIF(dce, mask_jad, aif, read_option, model_option) #run ETM 
    redcap_result.append({
            'field_name': 'MASK', 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value': 'JAD', 
        })
    for k in opt_params: 
        redcap_result.append({
            'field_name': model_option+'_'+k, 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value':opt_params[k], 
        })
except NameError: 
    print('JAD mask does not exist')

try: 
    opt_params = runtofts_popAIF(dce, mask_kwl, aif, read_option, model_option) #run ETM
    redcap_result.append({
            'field_name': 'MASK', 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value': 'KWL', 
        })
    for k in opt_params: 
        redcap_result.append({
            'field_name': model_option+'_'+k, 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value':opt_params[k], 
        })
except NameError: 
    print('KWL mask does not exist')

model_option = 'voxel'

try: 
    opt_params = runtofts_popAIF(dce, mask_jad, aif, read_option, model_option) #run ETM
    for k in opt_params: 
        redcap_result.append({
            'field_name': model_option+'_'+k, 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value':opt_params[k], 
        })
except NameError: 
    print()

try: 
    opt_params = runtofts_popAIF(dce, mask_kwl, aif, read_option, model_option) #run ETM
    for k in opt_params: 
        redcap_result.append({
            'field_name': model_option+'_'+k, 
            'record_id':description['PatientID'], 
            'baseline_arm_1': description['ReferringPhysician'],
            'value':opt_params[k], 
        })
except NameError: 
    print()

#save as json
with open(output+"/output.json", 'w') as outfile:
    outfile.write(json.dumps(redcap_result, indent=4)) 


