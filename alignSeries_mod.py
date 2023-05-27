from imagedata.series import Series
import numpy as np
#from showim import show3D_1, show3D_2
from scipy.interpolate import RegularGridInterpolator
#from read_b import *

def alignSeries_mod(dicomdirref, sortref, dicomdirmov, sortmov, movtype):

    # Default values
    #bvals = ''

    # -------------------- Load data ---------------------------------

    # Read the reference image
    imagepath = dicomdirref

    if len(sortref) == 0:
        imref = Series(imagepath)
    elif sortref == 'time':
        imref = Series(imagepath, sortref)
        timeline = imref.timeline
    else: 
        imref = Series(imagepath, sortref)

    dimref = imref.shape
 
    # Read the moving image that we want to resample
    imagepath = dicomdirmov

    if len(sortmov) == 0:
        immov = Series(imagepath)
    elif sortmov == 'time':
        immov = Series(imagepath, sortmov)
        timeline = immov.timeline
    else: 
        immov = Series(imagepath, sortmov)
    dimmov = immov.shape
   
    # Must to this in order to loop through the data
    if len(dimmov) == 3:
        immov = np.reshape(immov, [1, dimmov[0], dimmov[1], dimmov[2]])

    dimmov = immov.shape[1:]
    nch = immov.shape[0]
    
    if len(dimref) == 3: 
        # Make reference grid in voxel coordinates
        x = np.arange(0, dimref[0])
        y = np.arange(0, dimref[1])
        z = np.arange(0, dimref[2])
        [cx, cy, cz] = np.meshgrid(x, y, z, indexing='ij')
        nc = np.ones(np.prod(dimref)) #product 
        cref = np.asanyarray([cx.flatten(), cy.flatten(), cz.flatten(), nc], dtype='int')

        # Convert reference voxel coordinates to real coordinates
        xref = np.dot(imref.transformationMatrix, cref) #dot product

        # Generate voxel coordinates of reference image in the space of the moving image
        qinv = np.linalg.pinv(immov.transformationMatrix)
        cref2mov = np.dot(qinv, xref)

        # Only use the first three coordinates, the last is just ones. Interpolation method requires transpose()
        cref2mov = cref2mov[0:3, :].transpose()
    
        x = np.arange(0, dimmov[0])
        y = np.arange(0, dimmov[1])
        z = np.arange(0, dimmov[2])
        imreg = np.zeros([nch, dimref[0], dimref[1], dimref[2]], dtype='float')
        # Must convert to np array to slice the data
        immovnp = np.array(immov, dtype='float')
        for i in np.arange(nch):
            imh = immovnp[i,:,:,:]
            fnc = RegularGridInterpolator((x, y, z), imh, method='linear', bounds_error=False, fill_value=0)

            # Apply interpolator
            imh = fnc(cref2mov)
            imreg[i,:,:,:] = np.reshape(imh, dimref)

        # Squeeze
        imreg = np.squeeze(imreg)

    elif len(dimref) == 4: 

        x = np.arange(0, dimref[1])
        y = np.arange(0, dimref[2])
        z = np.arange(0, dimref[3])
        [cx, cy, cz] = np.meshgrid(x, y, z, indexing='ij')
        nc = np.ones(np.prod([dimref[1], dimref[2], dimref[3]]))
        cref = np.asanyarray([cx.flatten(), cy.flatten(), cz.flatten(), nc], dtype='int')                 

        # Convert reference voxel coordinates to real coordinates
        xref = np.dot(imref.transformationMatrix, cref)

        # Generate voxel coordinates of reference image in the space of the moving image
        qinv = np.linalg.pinv(immov.transformationMatrix)
        cref2mov = np.dot(qinv, xref)

        # Only use the first three coordinates, the last is just ones. Interpolation method requires transpose()
        cref2mov = cref2mov[0:3,:].transpose()

        x = np.arange(0, dimmov[0])
        y = np.arange(0, dimmov[1])
        z = np.arange(0, dimmov[2])

        imreg = np.zeros([1, dimref[1], dimref[2], dimref[3]], dtype='float')
        # Must convert to np array to slice the data
        immovnp = np.array(immov, dtype='float')
        imh = immovnp[0,:,:,:]
        fnc = RegularGridInterpolator((x, y, z), imh, method='linear', bounds_error=False, fill_value=0)

        # Apply interpolator
        imh = fnc(cref2mov)
        imreg = np.reshape(imh, dimref[1:4])

        # Squeeze
        imreg = np.squeeze(imreg)

    if movtype == 'mask':
        imreg = imreg > 0.5   

    imregout = Series(imreg)
    imregout = imregout.astype('float')

    return imref, imregout, dimref, dimmov, timeline