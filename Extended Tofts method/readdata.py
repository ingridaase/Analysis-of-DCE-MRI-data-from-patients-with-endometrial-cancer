from alignSeries_mod import * #The alignSeries function is now added to the imagedata library, see imagedata documentation https://imagedata.readthedocs.io/en/latest/UserGuide.html
from imagedata import Series
'''
This function reads the DICOM image and co-registrates to images with different grids.
'''

def readdata(mask_path, dce_path, option): 
    if option == 'mask2dce': 
        imref_dce, imreg_mask, dim_dce, dim_mask, timeline = alignSeries_mod(dce_path, 'time', mask_path, '', 'mask')
        MASK = imreg_mask
        DCE = imref_dce
    elif option == 'dce2mask': 
        imref_mask, imreg_dce, dim_mask, dim_dce, timeline = alignSeries_mod(mask_path, '', dce_path, 'time', '')
        DCE = imreg_dce
        MASK = imref_mask
    else: 
        print('NOT A VALID OPTION')

    return MASK, DCE, timeline
