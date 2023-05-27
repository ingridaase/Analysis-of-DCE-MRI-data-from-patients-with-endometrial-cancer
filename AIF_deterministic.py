import math
import pandas as pd
import numpy as np
import scipy.ndimage
from scipy.optimize import curve_fit
import skimage
import skimage.measure
import matplotlib.pyplot as plt
from imagedata.series import Series
import os 

'''
def show(signal, ax=None, show=False, label=None, color=None):
    if ax is None:
        fig, ax = plt.subplots(1)
        show = True
    ax.plot(signal, label=label, color=color)
    if label is not None:
        ax.legend(loc='upper right')
    if show:
        plt.show()
        plt.close()
'''

def relative_concentration_map(im, baseline, k, option):
    #print('This is computing {} signal concentration map'.format(
    #    option.upper())
    #)

    relim = np.zeros_like(im, dtype=np.float64)
    im0 = np.mean(im[:baseline], axis=0)
    for i in range(im.shape[0]):
        if option.upper().strip() == 'ABS':
            relim[i] = k * (im[i] - im0)
        elif option.upper().strip() == 'REL':
            relim[i] = k * (im[i] - im0) / im0
    relim[relim < 0] = 0
    return Series(relim, 'time', template=im, geometry=im)


#Step 1: Select 1% brightest voxels 
def select_brightest_voxels(dce, timesteps, percentile):
    assert timesteps>0 and timesteps<len(dce), "Timesteps {} out of range".format(timesteps)
    assert 0 < percentile < 100, "Percentile {} out of range".format(percentile)
    dce_sum = (np.sum(dce[:timesteps], axis=0, keepdims=True)/float(timesteps))[0]
    dce_pct = np.percentile(dce_sum, 100-percentile)
    brightest_voxels = dce_sum > dce_pct
    #dce_rgb = dce_sum.to_rgb()
    #dce_rgb[brightest_voxels] = (0xff, 0xff, 0)
    return brightest_voxels


#Step 2: Binary opening
def select_binary_opening(voxels):
    # morphological operation of a binary opening
    return scipy.ndimage.binary_opening(voxels)

#Step 3: Find timestep 
def find_most_peak_values(dce, mask, timesteps, ax1=None, ax2=None):
    ## Timestep of peak value for every voxel in mask
    peak_values = np.array(np.argmax(dce[:timesteps, mask], axis=0))
    #print('find_most_peak_values: peak_values type {}, nonzero {}, size {}, shape {}'.format(type(peak_values), np.count_nonzero(peak_values), peak_values.size, peak_values.shape))

    # Include voxels in first part of timeline only, do not include boundaries
    one_part = timesteps
    first_part_1 = np.bitwise_and(0 < peak_values, peak_values < one_part)
    
    first_part = np.zeros_like(dce[0], dtype=bool)
    first_part[mask] = first_part_1
    #print('find_most_peak_values: first_part', first_part.shape, first_part.sum())
    
    # Select most frequent peak timestep
    peak_value = np.argmax(dce[:timesteps, first_part], axis=0)
    hist, edges = np.histogram(peak_value, bins=one_part)
    max_timestep = np.argmax(hist)
    
    #print('find_most_peak_values: Max timestep: {}'.format(max_timestep))
    '''
    if ax1 is not None:
        #plt.subplot(221)
        ax1.plot(peak_values)
    if ax2 is not None:
        #plt.subplot(222)
        ax2.plot(hist)
        ax2.axvline(max_timestep, color='r')
    '''
    if max_timestep == 0 or max_timestep >= timesteps-1:
        # Do not include boundaries
        return None, None
    else:
        # Retain found timestep only
        retain_voxels = np.zeros_like(dce[0])
        retain_voxels[mask] = dce[max_timestep, mask]
        return max_timestep, retain_voxels

#Step 4: Fit Parker Population AIF
def find_promising_objects(regions, area_threshold):
    masks = []
    bbox = []
    list_of_index = []
    for num, x in enumerate(regions):
        area = x.area
        # convex_area = x.convex_area
        # if (num!=0 and (area>1000) and (convex_area/area <1.05) and (convex_area/area >0.95)):
        if num!=0 and (area>area_threshold):
            masks.append(regions[num].convex_image)
            bbox.append(regions[num].bbox)   
            list_of_index.append(num)
    return masks, bbox, list_of_index

def parker(x, N, t):
    # print('parker: x {} '.format(x.shape))
    A1, A2, T1, T2, sigma1, sigma2, alpha, beta, s, tau = x
    # print(A1,A2,T1,T2,sigma1,sigma2,alpha,beta,s,tau)

    term = np.zeros((3, N))
    term[0] = A1 / sigma1 / math.sqrt(2 * math.pi) * np.exp(-((t - T1) ** 2) / 2 / (sigma1 ** 2))
    term[1] = A2 / sigma2 / math.sqrt(2 * math.pi) * np.exp(-((t - T2) ** 2) / 2 / (sigma2 ** 2))
    term[2] = alpha * np.exp(-beta * t) / (1 + np.exp(-s * (t - tau)))

    # aif = np.sum(term, axis=0) * scaling
    aif = np.sum(term, axis=0)
    # return aif, term[0], term[1], alpha*np.exp(-beta*t), 1/(np.exp(-s*(t-tau)))
    return aif

def parker_cost_function(x, t, y):
    return parker(x, y.shape[0], t) - y

def estimate_parkers_model(measured, timeline):
    # timeline in minutes
    norm_factor = max(measured)
    y = measured / norm_factor

    A1 = 0.809
    A2 = 0.33
    T1 = 0.17046
    T2 = 0.365
    sigma1 = 0.0563
    sigma2 = 0.132
    alpha = 1.050
    beta = 0.1685
    s = 38.078
    tau = 0.483

    x0 = np.array((A1, A2, T1, T2, sigma1, sigma2, alpha, beta, s, tau))
    lb = np.array((0, 0, 0.1, 0.2, 1e-9, 1e-9, 1e-3, 0, 0, 0))
    ub = np.array((5, 5, 2, 2, 0.5, 0.7, 5.0, 1.17, 50, 1.5))
    xs = np.array((1, 1, 1, 1, 20, 10, 1, 1, 0.05, 1))
    result = scipy.optimize.least_squares(parker_cost_function, x0, bounds=(lb, ub), x_scale=xs,
                                    args=(timeline, y), verbose=1, ftol=1e-3, xtol=1e-3, gtol=1e-3)
    
    print('estimate_parkers_model: result.cost: {}'.format(result.cost))
    print('estimate_parkers_model: A1    : {:6.3f} {:6.3f}'.format(x0[0], result.x[0]))
    print('estimate_parkers_model: A2    : {:6.3f} {:6.3f}'.format(x0[1], result.x[1]))
    print('estimate_parkers_model: T1    : {:6.3f} {:6.3f}'.format(x0[2], result.x[2]))
    print('estimate_parkers_model: T2    : {:6.3f} {:6.3f}'.format(x0[3], result.x[3]))
    print('estimate_parkers_model: sigma1: {:6.3f} {:6.3f}'.format(x0[4], result.x[4]))
    print('estimate_parkers_model: sigma2: {:6.3f} {:6.3f}'.format(x0[5], result.x[5]))
    print('estimate_parkers_model: alpha : {:6.3f} {:6.3f}'.format(x0[6], result.x[6]))
    print('estimate_parkers_model: beta  : {:6.3f} {:6.3f}'.format(x0[7], result.x[7]))
    print('estimate_parkers_model: s     : {:6.3f} {:6.3f}'.format(x0[8], result.x[8]))
    print('estimate_parkers_model: tau   : {:6.3f} {:6.3f}'.format(x0[9], result.x[9]))

    return result
    
    
#Step 5: Erosion and dilation 

square = np.array([[1,1,1],
                   [1,1,1],
                   [1,1,1]])

def multi_dil(im, num, element=square):
    for i in range(num):
        im = skimage.morphology.dilation(im, element)
    return im
def multi_ero(im, num, element=square):
    for i in range(num):
        im = skimage.morphology.erosion(im, element)
    return im

#Step 6: Region growing 
def find_brightest_voxel(dce):
    assert dce.ndim == 3, "dce shall be 3-dim, is {}".format(dce.shape)
    _w = np.where(dce == np.max(dce))
    print('_w:', _w)
    try: 
        w0 = (int(_w[0]),)
        w1 = (int(_w[1]),)
        w2 = (int(_w[2]),)
    except TypeError:
        w0 = (int(_w[0][0]),)
        w1 = (int(_w[1][0]),)
        w2 = (int(_w[2][0]),)
    return w0, w1, w2

#Step 7: K-means clustering of time courses 
from sklearn.cluster import KMeans
def kmeans_cluster(df, n_clusters=1):
    """k-means cluster a given data frame into a specified number of
    clusters.
    It will then return a copy of the original data frame with those
    clusters appended in a column named Cluster.

    https://www.districtdatalabs.com/data-exploration-with-python-2
    """

    model = KMeans(n_clusters=n_clusters, random_state=1)
    clusters = model.fit_predict(df)
    #cluster_results = df.copy()
    #cluster_results['Cluster'] = clusters
    #return cluster_results
    return clusters

def summarize_clustering(results):
    """Count the number of objects that fall into each cluster and
    calculate the cluster means for each feature.
    It is going to merge the counts and means into a single data frame
    and then return that summary to us.

    https://www.districtdatalabs.com/data-exploration-with-python-2
    """

    cluster_size = results.groupby(['Cluster']).size().reset_index()
    cluster_size.columns = ['Cluster', 'Count']
    cluster_means = results.groupby(['Cluster'], as_index=False).mean()
    cluster_summary = pd.merge(cluster_size, cluster_means, on='Cluster')
    return cluster_summary

def find_deterministic_aif(dce): 

    experiment = {
        'dce': Series(dce, 'time'),
        'number_baseline': 5  # Length of baseline
    }
    
    patientID = experiment['dce'].patientID
    
    #Voxel volume 
    spacing = experiment['dce'].spacing
    z = spacing[0]
    y = spacing[1]
    x = spacing[2]
    pixel_size = x*y #mm^2
    voxel_size = x*y*z #mm^3
    ml = voxel_size*0.001 #milliliter
    
    # Show first volume of DCE
    #experiment['dce'][0].show(level=50, window=100)

    # Relative concentration
    k = 1
    dce_rel = relative_concentration_map(experiment['dce'], experiment['number_baseline'], k, 'abs')

    #Step 1: Find 1 percent brightest voxels during the 25 first timesteps
    brightest_voxels = select_brightest_voxels(dce_rel, timesteps=20, percentile=2)
    #print('Portion of brightest voxels: {:.1f}%'.format(100*brightest_voxels.sum() / dce_rel[0].size))
    #dce_rel[20,5].show(brightest_voxels[5])
    
    #Step 2: Binary Opening 
    binary_opening = select_binary_opening(brightest_voxels)
    #print(binary_opening.shape)
    #print('Portion of binary opening voxels: {:.1f}%'.format(100*binary_opening.sum() / dce_rel[0].size))
    #fig, ax = plt.subplots(1)
    #plt.imshow(binary_opening[5])

    #Step 3: Find timestep with most peak values
    #fig, axs = plt.subplots(2, 2)
    max_timestep, retain_voxels = find_most_peak_values(dce_rel, binary_opening, timesteps=20)
    if max_timestep is None:
        max_timestep, retain_voxels = find_most_peak_values(dce_rel, binary_opening, timesteps=30)
    #print('Portion of retained voxels 1: {:.1f}%'.format(100*np.count_nonzero(retain_voxels) / dce_rel[0].size))
    '''
    try: 
        axs[1,0].imshow(retain_voxels[5])
        plt.show()
    except TypeError:
        print('No retained voxels')
    '''
    
    #Step 4: Fit Parker population AIF
    #print('binary_opening', binary_opening.shape, binary_opening.dtype)
    label_im, num_labels = skimage.measure.label(binary_opening, return_num=True, connectivity=3)

    regions = skimage.measure.regionprops(label_im)
    '''
    print('label_im', label_im.shape, label_im.dtype, num_labels)
    print('Regions: {}'.format(len(regions)))
    fig, ax = plt.subplots(1, 1)
    ax.imshow(label_im[5])
    plt.show()
    plt.close()
    '''
    properties = ['area','convex_area','bbox_area', 'extent',
                  'mean_intensity', 'solidity', #'eccentricity',
                  #'orientation'
        ]
    df = pd.DataFrame(skimage.measure.regionprops_table(label_im, binary_opening, 
                 properties=properties))
    df['area'] = df['area'].apply(lambda x: x*ml)
    df.rename(columns = {'area':'area [ml]'}, inplace = True)
    #print(df)
    
    area_threshold = 1000
    #print(df[df['area [mm^3]'] > area_threshold])
    
    masks, bbox, list_of_index = find_promising_objects(regions, area_threshold)
    #print('masks:', masks)
    count = len(masks)
    painting = dce_rel[0].to_rgb()
    #print('Painting:', painting.shape, painting.dtype)
    try: 
        fig, ax = plt.subplots(3, int(count//3), figsize=(15,8))
    except ValueError: 
        fig, ax = plt.subplots(3, 1, figsize=(15,8))
        
    for axis, box, mask in zip(ax.flatten(), bbox, masks):
        red  =  painting[...,0][box[0]:box[3], box[1]:box[4], box[2]:box[5]] * mask
        green = painting[...,1][box[0]:box[3], box[1]:box[4], box[2]:box[5]] * mask
        blue  = painting[...,2][box[0]:box[3], box[1]:box[4], box[2]:box[5]] * mask
        mid = len(red)//2
        image = Series(np.zeros(red.shape+(3,), dtype=red.dtype))
        image[..., 0] = red
        image[..., 1] = green
        image[..., 2] = blue
        #axis.imshow(image[mid]*10)
    #plt.show()
    
    rgb_mask = np.zeros_like(label_im)
    for x in list_of_index:
        rgb_mask += (label_im==x+1).astype(int)
    red  =  painting[...,0] * rgb_mask
    green = painting[...,1] * rgb_mask
    blue  = painting[...,2] * rgb_mask
    image = dce_rel[0].to_rgb()
    #print('image', image.shape, image.dtype)
    image[..., 0] = red
    image[..., 1] = green
    image[..., 2] = blue
    #plt.imshow(image[5]*10)
    #plt.show()
    
    time_course = {}
    #fig, ax = plt.subplots(len(list_of_index), figsize=(15,8))
    for i, x in enumerate(list_of_index):
        rgb_mask = (label_im==x+1).astype(int)
        time_course[x] = np.sum(dce_rel, axis=(1,2,3), where=label_im==x+1)/np.count_nonzero(label_im==x+1)
        #show(time_course[x], ax=ax[i], label=x)
    #plt.show()
    
    # Fit Parker population AIF
    fit = {}
    df['cost'] = math.nan
    df['T1'] = math.nan
    df['T2'] = math.nan
    #fig, ax = plt.subplots(len(list_of_index), figsize=(15,8))
    for i, x in enumerate(list_of_index):
        #print('Fit time course {}:'.format(x))
        half = len(dce_rel.timeline)//2
        print(half)
        # solution = estimate_parkers_model(time_course[x], dce_rel.timeline/60)
        solution = estimate_parkers_model(time_course[x][:half], dce_rel.timeline[:half]/60)
        #print('solution', solution.cost, solution.success)
        fit[x] = solution
        #show(time_course[x], ax=ax[i], label=x)
        t_auc = np.sum(time_course[x][:half])
        p = parker(fit[x].x, half, fit[x].fun.timeline[:half]/60)
        p_auc = np.sum(p)
        #show(p*t_auc/p_auc,ax=ax[i], label=x)
        df.at[x, 'cost'] = solution.cost
        df.at[x, 'T1'] = solution.x[2]
        df.at[x, 'T2'] = solution.x[3]
    #plt.show()
    #plt.close()
    
    print(df[df['cost'] > 1e-3])
    print(df)
    i = df['cost'].idxmin()
    best_fit = i
    try: 
        new_df = df.iloc[[i]]
    except IndexError: 
        new_df = df
    print(new_df)
        
    
    #Step 5: Erosion and dialtion 
    mask = (label_im==best_fit+1).astype(int)
    #print(type(mask), mask.shape, mask.dtype, mask.max())

    #multi_dilated = multi_dil(mask, num=3, element=skimage.morphology.cube(3))
    #multi_eroded = multi_ero(multi_dilated, num=3, element=skimage.morphology.cube(3))
    multi_eroded = multi_ero(mask, num=1, element=skimage.morphology.cube(3))
    multi_dilated = multi_dil(multi_eroded, num=1, element=skimage.morphology.cube(3))
    #print('multi_eroded', multi_eroded.max())
    #print('multi_dilated', multi_dilated.max())
    #plt.subplot(221)
    #plt.imshow(mask[5]*100)
    #plt.subplot(222)
    #plt.imshow(multi_eroded[5]*100)
    #plt.subplot(223)
    #plt.imshow(multi_dilated[5]*100)
    # mask = multi_eroded
    #plt.tight_layout()
    #plt.show()
    
    tc = np.sum(dce_rel, axis=(1,2,3), where=multi_dilated==True)/np.count_nonzero(multi_dilated)
    print(tc)
    print(tc[:half])
    try: 
        solution = estimate_parkers_model(tc[:half], dce_rel.timeline[:half]/60)
        if solution.cost < df['cost'].min():
            #print('Selecting dilated mask')
            mask = multi_dilated
            cost_after_dilation = solution.cost
        else:
            #print('Keeping mask')
            cost_after_dilation = df['cost'].min()
    except ValueError: 
        #print('Keeping mask')
        cost_after_dilation = df['cost'].min()
    
    new_df['cost after dilation'] = cost_after_dilation
    
    #Step 6: Region growing 
    #print('mask', mask.shape, mask.dtype, mask.max(), np.count_nonzero(mask))
    _sum = np.sum(dce_rel[:10] * mask, axis=0)/10.
    seed_point = find_brightest_voxel(_sum)
    print('seed', seed_point)
    #print('_sum', _sum.shape, _sum.dtype)

    _sum_in_mask = np.extract(mask==1, _sum)
    # print('_sum_in_mask', _sum_in_mask.shape, np.count_nonzero(mask))

    brightest_intensity = float(_sum[seed_point])
    _quantile = np.quantile(_sum_in_mask, 0.8)
    _sum_in_quantile = np.extract(_sum_in_mask>_quantile, _sum_in_mask)
    #plt.subplots(1, 1)
    #plt.hist(_sum_in_quantile)
    #print('quantile', _quantile)
    _mean = np.mean(_sum_in_mask, where=_sum_in_mask>_quantile)
    _std = np.std(_sum_in_mask, where=_sum_in_mask>_quantile)
    threshold = _mean - _std
    #print('brightest', brightest_intensity, 'mean', _mean, 'std', _std, 'threshold', threshold)
    #plt.show()
    
    # Flood fill region growing from seed_point with threshold
    #print('grow: threshold', threshold)
    tolerance=brightest_intensity-threshold
    #print('grow: tolerance', tolerance)
    #print('grow: seed_point', seed_point)
    seed_point=(seed_point[0][0], seed_point[1][0], seed_point[2][0])
    #print('grow: seed_point', seed_point)
    m = skimage.morphology.flood(_sum, seed_point, tolerance=tolerance)
    #print('grow: m', m.shape, m.dtype, m.max(), np.count_nonzero(m))
    '''
    plt.subplots(1, 2)
    plt.subplot(121)
    plt.imshow(mask[5])
    plt.subplot(122)
    plt.imshow(m[5])
    '''
    tc = np.sum(dce_rel, axis=(1,2,3), where=m==True)/np.count_nonzero(m)
    solution = estimate_parkers_model(tc[:half], dce_rel.timeline[:half]/60)
    #print('Cost after region growing: {:.4f}, before region growing: {:.4f}'.format(solution.cost, cost_after_dilation))
    if solution.cost < cost_after_dilation:
        #print('Selecting region grow mask')
        mask = m
        cost_after_region_growing = solution.cost
    else:
        #print('Keeping mask')
        cost_after_region_growing = cost_after_dilation
    new_df['cost after region growing'] = cost_after_region_growing 
    
    #K-means clustering of time courses 
    #plt.subplots(1, 1)
    #plt.imshow(mask[5])
    timecourses = []
    for s in range(dce_rel.shape[1]):
        for r in range(dce_rel.shape[2]):
            for c in range(dce_rel.shape[3]):
                if mask[s,r,c]:
                    timecourses.append(dce_rel[:,s,r,c])
    #print(len(timecourses))
    #plt.show()
    #plt.close()
    
    cluster_results = kmeans_cluster(timecourses, 1)
    number_of_clusters = cluster_results.max()+1
    #print(len(cluster_results), number_of_clusters)
    hist, edges = np.histogram(cluster_results, range=(0,number_of_clusters), bins=number_of_clusters)
    #print(hist)
    #print(edges)
    
    tc = {}
    idx = {}
    for x in range(number_of_clusters):
        tc[x] = np.empty([hist[x],dce_rel.shape[0]], dtype=dce_rel.dtype)
        #print(x, tc[x].shape)
        idx[x] = 0
    i = 0
    cmask = np.zeros_like(dce_rel[0], dtype=int)
    for s in range(dce_rel.shape[1]):
        for r in range(dce_rel.shape[2]):
            for c in range(dce_rel.shape[3]):
                if mask[s,r,c]:
                    cluster = cluster_results[i]
                    ind = idx[cluster]
                    tc[cluster][ind,:] = dce_rel[:,s,r,c]
                    cmask[s,r,c] = cluster+1
                    idx[cluster] += 1
                    i +=1
    # Average clusters
    tcc = {}
    baseline = np.zeros(number_of_clusters)
    fullline = np.zeros(number_of_clusters)
    for x in range(number_of_clusters):
        assert idx[x] == hist[x], "idx {} disagree with hist {}".format(idx[x], hist[x])
        tcc[x] = np.sum(tc[x], axis=0)/tc[x].shape[0]
        b = experiment['number_baseline']
        baseline[x] = np.sum(tcc[x][:b])/b
        fullline[x] = np.sum(tcc[x])/len(tcc[x])
    best_fit = np.argmin(fullline)
    #print(best_fit)
    '''
    import matplotlib as mpl
    cc = ['r', 'g', 'b', 'c', 'y']
    cmap = mpl.colors.ListedColormap(cc)
    
    fig, ax = plt.subplots(figsize=(15,8))
    for x in range(number_of_clusters):
        if x == best_fit:
            show(tcc[x], ax=ax, label='{} *'.format(x), color=cc[x])
        else:
            show(tcc[x], ax=ax, label=x, color=cc[x])
    plt.show()
    plt.close()
    '''
    
    new_df['patientID'] = patientID
    cost_dict = new_df.to_dict()
    
    result_path = 'H:/data/Results/AutoAIF/2percent/AIF_plots/'
    
    result = best_fit
    Cb = tcc[best_fit]
    timeline = experiment['dce'].timeline
    plt.plot(timeline,Cb)
    plt.title('Arterial Input Function')
    plt.savefig(os.path.join(result_path, patientID))
    plt.clf() #Clear figure
    
    return Cb, cost_dict