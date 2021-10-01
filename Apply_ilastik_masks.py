#!/usr/bin/env python
# coding: utf-8

# ### Apply ilastik masks and save cell-level data with metadata
# 

# In[2]:


#setup libraries
import pandas as pd
import numpy as np
import os, re, glob, sys
from skimage import io, morphology, segmentation, measure
import h5py


# Instead of using Cellprofiler, convert the pixel masks to nuclei masks and apply them to the green images

# In[15]:


pipeline_name = 'PI' # python + ilastik
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
plateID = sys.argv[1]
well = sys.argv[2]
input_files_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
image_stack_paths = sorted(glob.glob(input_files_path+"/"+plateID+"_RGP_"+well+"*stack.tif"))
mask_paths = sorted(glob.glob(input_files_path+"/"+plateID+"_RGP_"+well+"*Segmentation.h5"))
mainpath = os.path.join(data_path,plateID,"Analysis",pipeline_name)


# In[16]:


ch1_name = 'NR'
ch2_name = 'CC'
cyto_expansion = 5
minimum_nuclear_radius = 2
minimum_nuclear_area = 3.14*minimum_nuclear_radius**2


# Get the green raw image files and apply the mask to extract the intensity data. Extract morphology values from the masks.

# In[19]:


results = []
for mask_path in mask_paths: #Each mask_path is a sequence of images
    print("processing files in "+mask_path)
    #load the pixel masks
    f = h5py.File(mask_path, 'r')
    mask_dataset = f['exported_data']
    masks_stack = np.stack(mask_dataset)
    nuclei_masks_raw = masks_stack == 4 # set nuclei pixels to True and the rest to False

    #load the corresponding green images
    field = re.findall("_[1-9]_",mask_path)[0]
    field_num = re.findall("[1-9]", field)[0]

    g_data_paths = glob.glob(os.path.join(data_path,plateID,well,"field_"+field_num,"*_G_*.tif"))
    img_g_ic = io.imread_collection(g_data_paths) # 3 dimensions : frames x width x height sorted by name
    img_gs = np.stack(img_g_ic)
    
    #read in the corresponding data file names to get time slices
    with open(re.sub("stack_Simple Segmentation.h5","filenames.txt", mask_path)) as f:
        filenames = f.readlines()
    
    for img_num, image in enumerate(nuclei_masks_raw[:,:,:,0]): #process each image in the mask_path sequence of images

        # open masks to delete small regions
        nuclei_masks_open = morphology.binary_opening(image, selem=morphology.disk(2))     

        # label the masks with unique integers starting at 0
        nuclei_masks_all = measure.label(nuclei_masks_open)
        #read in filenames to get time slice data
        #need to ensure there are nuclei pixels to process
        if np.amax(nuclei_masks_all) >0: #Only process if there is at least one mask
            nuclei_g = measure.regionprops_table(nuclei_masks_all, intensity_image = img_gs[img_num], properties=('label', 'area'))

            #remove masks too small to be a nucleus
            indices_to_keep = np.array([x if x-1 in np.argwhere(nuclei_g['area']>minimum_nuclear_area)
                                    else 0 for x in range(nuclei_g['label'].max()+1)])
            nuclei_masks = indices_to_keep[nuclei_masks_all]
            nuclei_g = measure.regionprops_table(nuclei_masks, intensity_image = img_gs[img_num], properties=('label', 'area','eccentricity',
                                                                                                                          'mean_intensity','max_intensity','min_intensity'))
            #expand the masks to get cytoplasmic regions
            nuclei_boundaries = segmentation.find_boundaries(nuclei_masks, mode='thick')*nuclei_masks
            nuclei_expansions = segmentation.expand_labels(nuclei_masks, cyto_expansion) - nuclei_masks + nuclei_boundaries
            nuclei_exp_g = measure.regionprops_table(nuclei_expansions, intensity_image = img_gs[img_num],
                                                     properties=('label','mean_intensity','max_intensity','min_intensity'))

            # turn results into a dataframe
            nuclei_g_data = pd.DataFrame(nuclei_g)
            nuclei_g_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_g_data.columns if col not in ['label']}, inplace=True)

            nuclei_exp_g_data = pd.DataFrame(nuclei_exp_g)
            nuclei_exp_g_data.rename(columns={col: 'Cyto_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_exp_g_data.columns if col not in ['label']}, inplace=True)

            # add an image number and collect the data                                                                                                             
            nuclei_g_data['image'] = img_num+1

            #Calculate ratio of ch2 cyto to nuclei intensities
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_mean_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_mean_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_mean_intensity']
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_max_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_max_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_max_intensity']
            nuclei_exp_g_data['Cell_'+pipeline_name+'_' +ch2_name+'_min_intensity_ratio'] = nuclei_exp_g_data['Cyto_'+pipeline_name+'_' +ch2_name+'_min_intensity']/nuclei_g_data['Nuclei_'+pipeline_name+'_' +ch2_name+'_min_intensity']

            # add the well and field values to the dataframe
            nuclei_g_data['well'] = well
            nuclei_g_data['field'] = field_num
            #add the time slice to the dataframe
            nuclei_g_data['time_slice'] = re.findall("[0-9d]{3}[0-9h]{3}[0-9m]{3}", filenames[img_num])[0]

            #concatenate the dataframes
            df_all = pd.concat([nuclei_g_data, nuclei_exp_g_data], axis=1, join="outer")
            #append this image's data to the rest of the data
            results.append(df_all)
    #Done processing all of the images in the sequence so concatenate all of the dataframes
    results_pd = pd.concat(results)

#concatenate all of the results from the sequences in the well
l0 = pd.concat(results)

    #Save mask image
    #mask_filename = data_path+plateID+'/'+well+'/'+field+'/output_stacks/'+well+'_'+field+'_image'+str(img_num)+'_nuclei_masks.tif'
        #io.imsave(mask_filename, nuclei_masks.astype('uint16'))
        #cyto_mask_filename = data_path+plateID+'/'+well+'/'+field+'/output_stacks/'+well+'_'+field+'_image'+str(img_num)+'_cyto_masks.tif'
        #io.imsave(cyto_mask_filename, nuclei_expansions.astype('uint16'))


# If the metadata files exists, combine the data with the experimental metadata and write out a level 1 file. Otherwise, write out a level 0 file so we don't waste the results.

# In[24]:


#If the metadata exists, join it to the data and write out as a level 1 file
metadata_filename = os.path.join(data_path,plateID,"metadata",plateID+".xlsx")
if not os.path.exists(mainpath):
    os.makedirs(mainpath, exist_ok=True)
if os.path.exists(metadata_filename):
    print("adding "+plateID+" metadata")
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
    
    #merge data and metadata on well values
    l1= pd.merge(l0, metadata, how="left", on=["well"])
    l1.to_csv(os.path.join(mainpath,plateID+'_'+well+'_level_1.csv'))
else:
    print("no metadata file for "+plateID+" so creating level 0 file")
    l0.to_csv(os.path.join(mainpath,plateID+'_'+well+'_level_0.csv'))


# In[ ]:




