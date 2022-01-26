#!/usr/bin/env python
# coding: utf-8

# In[1]:


#setup libraries
import numpy as np
import pandas as pd
import os, re, glob, sys
from skimage import io, filters, util, segmentation, morphology, measure, restoration, exposure
import matplotlib.pyplot as plt
#import imageio
from cellpose import models, plot
from pystackreg import StackReg
from scipy import stats
import tifffile


# In[2]:


pipeline_name = "PC" #python and cellpose
ch1_name = 'NR'
ch2_name = 'CC'
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
#plateID = 'AU02001'
#sys.argv[1] = /home/exacloud/gscratch/HeiserLab/images/AU02001/Analysis/PC/intermediate_files/AU02001_R_A1_1_reg_stack.tif
plateID = sys.argv[1]
#well_index = 1
well_index = int(sys.argv[2]) - 1

output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files/")
transformation_path = os.path.join(output_path,"transformations")
#Only use subdirectories within the selected well
well_directory = sorted(glob.glob(os.path.join(data_path+plateID,"[A-Z][1-9]")))[well_index]
well = re.findall("[A-Z][1-9]$",well_directory)[0]
subdirectories = sorted(glob.glob(os.path.join(well_directory,"field_[1-9]")))


# Only create a job if there is not a registered red channel stack

# In[3]:


flourescent_scaler = 255/4095 #rescale from 12 to 8 bits

for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    # Only process the field-level image files if there is no registered stack
    if not os.path.exists(reg_filename):
        print("registering R stack in "+subdir)
        #load and prepare red, green and phase channels. Scale for 8 bits but these are uint16 data types
        r_data_paths = glob.glob(os.path.join(subdir,"*_R_*m.tif"))
        r_time_slices = set()
        for data_paths in r_data_paths:
            r_time_slices.add(re.findall("..d..h..m", data_paths)[0])
        g_data_paths = glob.glob(os.path.join(subdir,"*_G_*m.tif"))
        g_time_slices = set()
        for data_paths in g_data_paths:
            g_time_slices.add(re.findall("..d..h..m", data_paths)[0])
        p_data_paths = glob.glob(os.path.join(subdir,"*_P_*m.tif"))
        p_time_slices = set()
        for data_paths in p_data_paths:
            p_time_slices.add(re.findall("..d..h..m", data_paths)[0])
        complete_time_slices = r_time_slices & g_time_slices & p_time_slices
        r_data_paths_c = []
        g_data_paths_c = []
        p_data_paths_c = []
        for time_slice in complete_time_slices:
            r_data_paths_c.append(os.path.join(data_path+plateID,well,field,plateID+"_R_"+well+"_"+field_num+"_"+time_slice+".tif"))
            g_data_paths_c.append(os.path.join(data_path+plateID,well,field,plateID+"_G_"+well+"_"+field_num+"_"+time_slice+".tif"))
            p_data_paths_c.append(os.path.join(data_path+plateID,well,field,plateID+"_P_"+well+"_"+field_num+"_"+time_slice+".tif"))
        img_r_ic = io.imread_collection(r_data_paths_c) # 3 dimensions : frames x width x height
        img_rs = np.stack(img_r_ic)*flourescent_scaler

        img_g_ic = io.imread_collection(g_data_paths_c) # 3 dimensions : frames x width x height
        img_gs = np.stack(img_g_ic)*flourescent_scaler

        img_p_ic = io.imread_collection(p_data_paths_c) # 3 dimensions : frames x width x height
        img_ps = np.stack(img_p_ic)
        
        #register the R stack using transformation
        sr = StackReg(StackReg.TRANSLATION)
        # register each frame to the previous (already registered) one
        tmats = sr.register_stack(img_rs, reference='previous')
        if not os.path.exists(transformation_path):
            os.makedirs(transformation_path)
        np.save(os.path.join(transformation_path,plateID+"_"+well+"_"+field+"_transformation_matrices.npy"), tmats)
        
        # transform stack using the tmats loaded from file
        img_rs_reg = sr.transform_stack(img_rs, tmats = tmats)
        img_gs_reg = sr.transform_stack(img_gs, tmats = tmats)
        img_ps_reg = sr.transform_stack(img_ps, tmats = tmats)

        img_c = np.stack([img_rs.astype('B'), img_gs.astype('B'), img_ps.astype('B')], axis = -1) 
    
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        io.imsave(reg_filename, img_rs_reg, plugin='tifffile')
        io.imsave(reg_filename.replace("_R_","_G_"), img_gs_reg, plugin='tifffile')
        io.imsave(reg_filename.replace("_R_","_P_"), img_ps_reg, plugin='tifffile')


# segment each registered stack using cellpose
#     

# In[4]:


#If the metadata exists, join it to the data and write out as a level 1 file
metadata_filename = os.path.join(data_path,plateID,"metadata",plateID+".xlsx")

if os.path.exists(metadata_filename):
    md_all = pd.read_excel(metadata_filename, engine='openpyxl', dtype={'Drug1Concentration': str, 'Drug2Concentration': str})
    
    #remove unwanted columns read in from the excel files
    r = re.compile("Unnamed.*")
    columns_to_drop = list(filter(r.match, md_all.columns)) 
    metadata = md_all.drop(columns = columns_to_drop)
    
    #match metadata and data well labels format
    metadata['column'] = [re.sub(r'[0-9]*', '', Well) for Well in metadata['Well']]
    metadata['row'] = [re.sub(r'[A-Z]', '', Well) for Well in metadata['Well']]
    metadata['row'] = [re.sub(r'\A0', '', row) for row in metadata['row']]
    metadata['well'] = metadata['column'] + metadata['row']
    


# In[ ]:


diameter = 13.2
flow_threshold =1.0
mask_threshold=0
min_size=75
resample = True
cyto_expansion = 5
minutes_between_images = 30

#call cellpose on each image to segment the nuclei
# DEFINE CELLPOSE MODEL
model = models.Cellpose(gpu=True, model_type='cyto')

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# will use channel R = 1 as nuclear channel only
channels = [0,1]

for subdir in subdirectories:
    results = [] #collect results for one field in the well
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    l0_filename = os.path.join(output_path,plateID+"_"+well+"_"+field+"_level_0.csv")
    if not os.path.exists(l0_filename): #Only segment if no level 0 file
        img_rs_reg = io.imread(reg_filename)
        img_gs_reg = io.imread(reg_filename.replace("_R_","_G_"))
        img_ps_reg = io.imread(reg_filename.replace("_R_","_P_"))
                         
        mask_images = []
        for i, image in enumerate(img_rs_reg):
            print("processing "+reg_filename+" index "+str(i))
            #Segment using the nuclear signal in the red channel
            background = restoration.rolling_ball(image, radius=50)
            image_bc = exposure.adjust_log(image-background, gain=1, inv=False)
            image_rf = filters.median(image_bc, morphology.disk(1))

            # create masks with cellpose 
            masks, flows, styles, diams = model.eval(image_rf,
                                             diameter=diameter,
                                             flow_threshold=flow_threshold,
                                             mask_threshold=mask_threshold,
                                             channels=channels,
                                             min_size=min_size,
                                             resample = resample)
            #measure reporter intensity and nuclear morphology, texture
            nuclei = measure.regionprops_table(measure.label(masks), intensity_image=image,
                                               properties=('label', 'area', 'bbox_area', 'convex_area', 'centroid', 'eccentricity',
                                                           'equivalent_diameter', 'extent', 'feret_diameter_max','filled_area',
                                                           'major_axis_length', 'minor_axis_length', 'moments_hu', 'perimeter',
                                                           'perimeter_crofton', 'solidity', 'mean_intensity', 'max_intensity',
                                                           'min_intensity'))
            #Save mask image
            mask_images.append(masks)
    
            #expand the masks to get cytoplasmic regions
            nuclei_boundaries = segmentation.find_boundaries(masks, mode='thick')*masks
            nuclei_expansions = segmentation.expand_labels(masks, cyto_expansion) - masks + nuclei_boundaries
    
            # measure nuclear and cytoplasmic intensities and textures in the green channel
            nuclei_g = measure.regionprops_table(measure.label(masks), intensity_image=img_gs_reg[i],
                                               properties=('label',
                                                            'mean_intensity','max_intensity','min_intensity'))
            nuclei_exp_g = measure.regionprops_table(measure.label(nuclei_expansions), intensity_image=img_gs_reg[i],
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
            well = re.findall('_[A-Z][0-9]+_',reg_filename)[0]
            well = re.sub('_','', well)
            nuclei_data['well'] = well
            field = re.findall('_[0-9]+_',reg_filename)[0]
            field = re.sub('_','', field)
            nuclei_data['field'] = field
            nuclei_data['slice'] = i
            elapsed_minutes = i*minutes_between_images #assumes time slice numbering starts at 1
            day = np.floor(elapsed_minutes/(24*60)).astype(int)
            hour = np.floor((elapsed_minutes-day*(24*60))/60).astype(int)
            minute = np.floor(elapsed_minutes-day*(24*60)-hour*60).astype(int)
            day = str(day).zfill(2)
            hour = str(hour).zfill(2)
            minute = str(minute).zfill(2)
            nuclei_data['time_slice'] = day+"d"+hour+"h"+minute+"m"

            #Calculate ratio of ch2 cyto to nuclei intensities
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_mean_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_mean_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_mean_intensity']
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_max_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_max_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_max_intensity']
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_min_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_min_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_min_intensity']

            #concatenate the dataframes from the different channels
            df_all = pd.concat([nuclei_data, nuclei_g_data, nuclei_exp_g_data], axis=1, join="outer")
            # append the dataframe to the results list
            results.append(df_all)
            
            #Save mask image
            mask_filename = reg_filename.replace("_reg_stack.tif", "_masks_stack.png")
            tifffile.imwrite(mask_filename, np.array(mask_images), imagej=True)
    
    #concatenate all of the results from the sequences in the well
    l0 = pd.concat(results)
    
    if os.path.exists(metadata_filename):
        #merge data and metadata on well values
        l1= pd.merge(l0, metadata, how="left", on=["well"])
        l1.to_csv(l0_filename.replace('level_0','level_1'))
    else:
        print("no metadata file for "+plateID+" so creating level 0 file")
        l0.to_csv(l0_filename)
            


# In[ ]:




