from alignSeries_mod import *
from imagedata.series import Series

#Fjernet mulitheten for å registrere fra dce til vibe, usikker på om jeg trenger denne da det tar lenger tid 
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