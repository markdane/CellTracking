#!/usr/bin/env python
# coding: utf-8

# ### AU565 Live Cell Images Segmentation
# 
# This notebook intiates processing of AU565  movies. The experiment was run in plate AU00901-01101 and images on the Incucyte microscope. Phase contast and two fluorescent channels were recorded.

#setup libraries
import os, re, sys, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, filters, segmentation, morphology, measure
from cellpose import models, plot
from scipy import stats

use_GPU = models.use_gpu()
print('>>> GPU activated? %d'%use_GPU)

pipeline_name = 'CL'
ch1_name = 'NR'
ch2_name = 'CC'
cyto_expansion = 5

data_path = '/home/exacloud/gscratch/HeiserLab/images/'
plateID = sys.argv[1]
well = sys.argv[2]
data_paths =glob.glob(os.path.join(data_path+plateID+'/'+well+'/','field_*/',"*_R_*m.tif"))
df = pd.DataFrame(data_paths, columns =['path']).sort_values(by='path', ignore_index = True)

#call cellpose on each image to segment the nuclei
# DEFINE CELLPOSE MODEL
# model_type='cyto' or model_type='nuclei'
model = models.Cellpose(gpu=use_GPU, model_type='cyto2')

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
channels = [0,1]
    
results = []
#Segment using the nuclear signal in the grayscale image from the red channel
flow_threshold = .4
diameter = 16
cellprob_threshold=0.0
min_size=125
resample = True
# create an index that allows us to go through each pair of images
for path_name in df['path']:
    path_name_p = str.replace(path_name, '_R_', '_P_')
    path_name_g = str.replace(path_name, '_R_', '_G_')
    
    # load the image of the red nuclear channel and run a smoothing filter
    image_r = io.imread(path_name)
    image_rf = filters.median(image_r, selem=morphology.disk(2))
    
    # load the phase image for display purposes
    image_p = io.imread(path_name_p)

    # load the phase image for display purposes
    image_g = io.imread(path_name_g)

    #Combine phase and nuclear images with nuceli in R 1 channel and phase in G 2 channel
    #image = np.dstack([image_rf, image_p])

    # create masks from filtered red channel with cellpose 
    masks, flows, styles, diams = model.eval(image_rf, diameter=diameter, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold, channels=channels, min_size=min_size, resample = resample)

    #measure reporter intensity and nuclear morphology, texture
    nuclei = measure.regionprops_table(measure.label(masks), intensity_image=image_r,
                                           properties=('label',
                                                       'area','bbox_area','convex_area','centroid','eccentricity','equivalent_diameter','extent','feret_diameter_max','filled_area',
                                                        'major_axis_length','minor_axis_length','moments_hu','perimeter','perimeter_crofton','solidity',
                                                        'mean_intensity','max_intensity','min_intensity'))
    #Save mask image
    mask_filename = path_name.replace(".tif", "_masks.png")
    io.imsave(mask_filename, masks)
    
    #expand the masks to get cytoplasmic regions
    nuclei_boundaries = segmentation.find_boundaries(masks, mode='thick')*masks
    nuclei_expansions = segmentation.expand_labels(masks, cyto_expansion) - masks + nuclei_boundaries
    
    # measure nuclear and cytoplasmic intensities and textures in the green channel
    nuclei_g = measure.regionprops_table(measure.label(masks), intensity_image=image_g,
                                           properties=('label',
                                                        'mean_intensity','max_intensity','min_intensity'))
    nuclei_exp_g = measure.regionprops_table(measure.label(nuclei_expansions), intensity_image=image_g,
                                                 properties=('label',
                                                             'mean_intensity','max_intensity','min_intensity'))
    
    # turn results into a dataframe
    nuclei_data = pd.DataFrame(nuclei)
    nuclei_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch1_name+'_'+col  for col in nuclei_data.columns if col not in ['label']}, inplace=True)
   
    nuclei_g_data = pd.DataFrame(nuclei_g)
    nuclei_g_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_g_data.columns if col not in ['label']}, inplace=True)
    

    nuclei_exp_g_data = pd.DataFrame(nuclei_exp_g)
    nuclei_exp_g_data.rename(columns={col: 'Cyto_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_exp_g_data.columns if col not in ['label']}, inplace=True)
       
    # recover the well and field values and add them to the dataframe
    well = re.findall('_[A-Z][0-9]+_',path_name)[0]
    well = re.sub('_','', well)
    nuclei_data['well'] = well
    field = re.findall('_[0-9]+_',path_name)[0]
    field = re.sub('_','', field)
    nuclei_data['field'] = field
    time_slice = re.findall('[a-z0-9]*.tif',path_name)[0]
    time_slice = re.sub('.tif','', time_slice)
    nuclei_data['time'] = time_slice
    
 #Calculate ratio of ch2 cyto to nuclei intensities
    nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_mean_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_'+ch2_name+'_mean_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_'+ch2_name+'_mean_intensity']
    nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_max_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_'+ch2_name+'_max_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_'+ch2_name+'_max_intensity']
    nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_min_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_'+ch2_name+'_min_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_'+ch2_name+'_min_intensity']
    
    #concatenate the dataframes from the different channels
    df_all = pd.concat([nuclei_data, nuclei_g_data, nuclei_exp_g_data], axis=1, join="outer")
    # append the dataframe to the results list
    results.append(df_all)
    print("processing "+path_name)
    
#concatenate all of the results
all_results = pd.concat(results)

mainpath = '../../../images/'+plateID+'/Analysis/CL/'
if not os.path.exists(mainpath):
    os.makedirs(mainpath, exist_ok=True)
all_results.to_csv(mainpath+plateID+'_'+well+'_level_0.csv')

