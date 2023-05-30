import numpy as np
import os
import matplotlib.pyplot as plt
from imagedata import Series
from skimage.morphology import erosion, dilation

#Functions 
from tofts_method import *
from readdata import *
from AIF_deterministic import *

'''
This script contains three functions for running ETM using three different methods for estimating the arterial input function. 
'''

#Run ETM using a manually annotated AIF.
def runtofts_manualAIF(dce, mask, AIF, read_option = 'mask2dce', model_option = 'average'): 

    #Read data 
    MASK, DCE, timeline = readdata(mask, dce, read_option)

    #Find voxel size
    spacing = DCE.spacing
    z = spacing[0]
    y = spacing[1]
    x = spacing[2]
    voxel_size = x*y*z #mm^3
    voxel_size = voxel_size*0.001 #ml
    
    #Subtract pre-contrast signal from signal
    DCE = DCE.astype('float32')
    S0 = np.mean(DCE[0:5,:,:,:], axis=0)
    dynim = DCE.copy()
    for k in range(len(timeline)):
        dynim[k,:,:,:] = DCE[k,:,:,:] - S0


    #Read AIF
    AIF_np = np.asarray(AIF)
    #Find mean values inside AIF mask 
    non_zero_in_AIF = np.nonzero(AIF_np) #index of nonzero elements in mask image
    AIF_vals = dynim[:, non_zero_in_AIF[1], non_zero_in_AIF[2], non_zero_in_AIF[3]]
    AIF_vals = np.asarray(np.mean(AIF_vals, axis=1), dtype='float')

    #Calculate plasma concentration from the value of hematocrit 
    cp = AIF_vals/(1-0.42) #42 percent hematocrit 

    if model_option == 'average': 
        #Finding mean concentration inside mask 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        in_mask = np.asarray(in_mask_vals, dtype='float')
        C = np.mean(in_mask, axis=1)
        tumor_vol = in_mask.shape[1]*voxel_size
        opt_params, C_tofts = run_tofts(timeline, C, cp)
        opt_params['tumor_volume'] = tumor_vol
    elif model_option == 'voxel': 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        in_mask_vals = np.asarray(in_mask_vals, dtype='float')
        tumor_vol = in_mask_vals.shape[1]*voxel_size
        C_tofts = np.zeros(in_mask_vals.shape)
        C = np.zeros(in_mask_vals.shape)
        params = []
        for i in range(in_mask_vals.shape[1]): 
            Cn = in_mask_vals[:,i]
            C[:,i] = Cn
            try: 
                opt_params, C_tofts_voxel = run_tofts(timeline, Cn, cp)
            #Handle possible Runtimeerror 
            except RuntimeError: 
                opt_params = {'ktrans':np.nan,
                  've': np.nan,
                  'vp':np.nan,
                  'kep':np.nan}
                C_tofts_voxel = np.nan
            params.append(opt_params)
            C_tofts[:,i] = C_tofts_voxel
        C_tofts = np.nanmean(C_tofts, axis=1)
        C = np.nanmean(C, axis=1)
        mean_dict = {}
        for key in params[0].keys():
            mean_dict[key] = np.nanmean([d[key] for d in params], axis=0)
        opt_params = mean_dict 
        opt_params['tumor_volume'] = tumor_vol

    return opt_params

#Run ETM using a automatic AIF.
def runtofts_autoAIF(dce_path, mask_path, read_option = 'mask2dce', model_option = 'average'): 
    
    #Read data 
    MASK, DCE, timeline = readdata(mask_path, dce_path, read_option)
    
    #Find pixel and voxel size
    spacing = DCE.spacing
    z = spacing[0]
    y = spacing[1]
    x = spacing[2]
    pixel_size = x*y
    voxel_size = x*y*z
    voxel_size = voxel_size*0.001 #ml
    
    Cb, Cb_params = find_deterministic_aif(DCE)
    
    #Subtract pre-contrast signal from signal
    DCE = DCE.astype('float32')
    S0 = np.mean(DCE[0:5,:,:,:], axis=0)
    dynim = DCE.copy()
    for k in range(len(timeline)):
        dynim[k,:,:,:] = DCE[k,:,:,:] - S0

    
    #Calculate plasma concentration from the value of red blood cells 
    cp = Cb/(1-0.42) #42 percent hematocrit (red blood cells)

    if model_option == 'average': 
        #Finding mean concentration inside mask 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        C = np.asarray(np.mean(in_mask_vals, axis=1), dtype='float')
        opt_params, C_tofts = run_tofts(timeline, C, cp)
        tumor_vol = in_mask_vals.shape[1]*voxel_size
        opt_params['tumor_volume'] = tumor_vol
    elif model_option == 'voxel': 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        in_mask_vals = np.asarray(in_mask_vals, dtype='float')
        tumor_vol = in_mask_vals.shape[1]*voxel_size
        C_tofts = np.zeros(in_mask_vals.shape)
        C = np.zeros(in_mask_vals.shape)
        params = []
        for i in range(in_mask_vals.shape[1]): 
            Cn = in_mask_vals[:,i]
            C[:,i] = Cn
            try: 
                opt_params, C_tofts_voxel = run_tofts(timeline, Cn, cp)
            except RuntimeError: 
                opt_params = {'ktrans': np.nan, 
                              've':np.nan, 
                              'vp':np.nan, 
                              'kep':np.nan}
                C_tofts_voxel = np.nan
            params.append(opt_params)
            C_tofts[:,i] = C_tofts_voxel
        C_tofts = np.mean(C_tofts, axis=1)
        C = np.mean(C, axis=1)
        mean_dict = {}
        for key in params[0].keys():
            mean_dict[key] = np.nanmean([d[key] for d in params], axis=0)
        opt_params = mean_dict 
        opt_params['tumor_volume'] = tumor_vol
    return opt_params

#Run ETM using a population-based AIF.
def runtofts_popAIF(dce_path, mask_path, aif, read_option = 'mask2dce', model_option = 'average'): 

    '''
    #Eoded or dilation tumor masks. 
    mask = Series(mask_path)
    eroded_mask = erosion(mask) 
    dilated_mask = dilation(mask)
    '''
    
    #Read data 
    MASK, DCE, timeline = readdata(mask_path, dce_path, read_option)
    
    patientID = DCE.patientID
    
    #Find pixel and voxel size
    spacing = DCE.spacing
    z = spacing[0]
    y = spacing[1]
    x = spacing[2]
    pixel_size = x*y
    voxel_size = x*y*z
    voxel_size = voxel_size*0.001 #ml
    
    cp = aif/(1-0.42)
    
    #Subtract pre-contrast signal from signal
    DCE = DCE.astype('float32')
    S0 = np.mean(DCE[0:5,:,:,:], axis=0)
    dynim = DCE.copy()
    for k in range(len(timeline)):
        dynim[k,:,:,:] = DCE[k,:,:,:] - S0

    if model_option == 'average': 
        #Finding mean concentration inside mask 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        C = np.asarray(np.mean(in_mask_vals, axis=1), dtype='float')
        opt_params, C_tofts = run_tofts(timeline, C, cp)
        tumor_vol = in_mask_vals.shape[1]*voxel_size
        opt_params['tumor_volume'] = tumor_vol
    elif model_option == 'voxel': 
        non_zero_in_mask = np.nonzero(MASK) #index of nonzero elements in mask image
        in_mask_vals = dynim[:, non_zero_in_mask[0], non_zero_in_mask[1], non_zero_in_mask[2]]
        in_mask_vals = np.asarray(in_mask_vals, dtype='float')
        tumor_vol = in_mask_vals.shape[1]*voxel_size
        C_tofts = np.zeros(in_mask_vals.shape)
        C = np.zeros(in_mask_vals.shape)
        params = []
        for i in range(in_mask_vals.shape[1]): 
            Cn = in_mask_vals[:,i]
            C[:,i] = Cn
            try: 
                opt_params, C_tofts_voxel = run_tofts(timeline, Cn, cp)
            except RuntimeError: 
                opt_params = {'ktrans': np.nan, 
                              've':np.nan, 
                              'vp':np.nan, 
                              'kep':np.nan}
                C_tofts_voxel = np.nan
            params.append(opt_params)
            C_tofts[:,i] = C_tofts_voxel
        C_tofts = np.mean(C_tofts, axis=1)
        C = np.mean(C, axis=1)
        mean_dict = {}
        for key in params[0].keys():
            mean_dict[key] = np.nanmean([d[key] for d in params], axis=0)
        opt_params = mean_dict 
        opt_params['tumor_volume'] = tumor_vol

    return opt_params
