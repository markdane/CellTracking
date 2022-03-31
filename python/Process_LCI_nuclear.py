#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
from PIL import Image, ImageSequence
import tifffile


# In[3]:


pipeline_name = "CKn" #cellpose and KIT tracking
ch1_name = 'NR'
ch2_name = 'CC'
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
#data_path = '/Users/dane/Documents/CellTrackingProjects/AU565/images/'
plateID = 'AU02001'
#plateID = sys.argv[1]
well_index = 21
#well_index = int(sys.argv[2])

output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files/")
transformation_path = os.path.join(output_path,"transformations")
tracking_path = os.path.join(output_path,'tracking/')

#Only use subdirectories within the selected well
well_directory = sorted(glob.glob(os.path.join(data_path+plateID,"[A-Z][1-9]")))[well_index-1]
well = re.findall("[A-Z][1-9]$",well_directory)[0]
subdirectories = sorted(glob.glob(os.path.join(well_directory,"field_[1-9]")))


# #### Register the image stacks
# If there is a registered red channel stack skip this step, otherwise:  
# Load the red, green and phase images  
# Delete images from any time slice that does not have a complete set of images  
# Rescale the fluorescent images from 12 to 8 bits  
# Calculate the transformations needed to register the red stack, correcting the translation only  
# Store the registration transformations  
# Use the transformations to register all three stacks  
# Save the registered stacks as 16 bit images  
# 
# TODO: See if it's possible to not rescale the images  

# In[ ]:


#fluorescent_scaler = 255/4095 #rescale from 12 to 8 bits
fluorescent_scaler = 1

for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    # Only process the field-level image files if there is no registered stack
    if not os.path.exists(reg_filename):
        print("Gathering image stacks for "+subdir)
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
        img_rs = np.stack(img_r_ic)*fluorescent_scaler

        img_g_ic = io.imread_collection(g_data_paths_c) # 3 dimensions : frames x width x height
        img_gs = np.stack(img_g_ic)*fluorescent_scaler

        img_p_ic = io.imread_collection(p_data_paths_c) # 3 dimensions : frames x width x height
        img_ps = np.stack(img_p_ic)
        
        #register the R stack using transformation
        print("Determining registration transformations for "+subdir)
        sr = StackReg(StackReg.TRANSLATION)
        # register each frame to the previous (already registered) one
        tmats = sr.register_stack(img_rs, reference='previous', axis = 0)
        if not os.path.exists(transformation_path):
            os.makedirs(transformation_path)
        np.save(os.path.join(transformation_path,plateID+"_"+well+"_"+field+"_transformation_matrices.npy"), tmats)
        
        # transform stack using the tmats loaded from file
        print("Using transformations to register all stacks for "+subdir)
        img_rs_reg = sr.transform_stack(img_rs, tmats = tmats)
        img_gs_reg = sr.transform_stack(img_gs, tmats = tmats)
        img_ps_reg = sr.transform_stack(img_ps, tmats = tmats)

        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        #crop each stack to the area that is common in all images after registration
        #use the min and max values in the transformation stack to define the common active area 
        x_max = tmats[:,0,2].max().astype("int")
        x_min = tmats[:,0,2].min().astype("int")
        y_max = tmats[:,1,2].max().astype("int")
        y_min = tmats[:,1,2].min().astype("int")
        
        #assume transformations are in pixels and crop images to exclude areas that are outside of any registered image
        x_axis_length = img_rs_reg.shape[2]
        y_axis_length = img_rs_reg.shape[1]
        print("Cropping to active areas for "+subdir)
        
        img_rs_reg_crop = img_rs_reg[:,-y_min:(y_axis_length-y_max),-x_min:(x_axis_length-x_max)]
        img_gs_reg_crop = img_gs_reg[:,-y_min:(y_axis_length-y_max),-x_min:(x_axis_length-x_max)]
        img_ps_reg_crop = img_ps_reg[:,-y_min:(y_axis_length-y_max),-x_min:(x_axis_length-x_max)]
        print("Saving stacks to disk for "+subdir)
        io.imsave(reg_filename, img_rs_reg_crop.astype(np.int16), plugin='tifffile', check_contrast=False)
        io.imsave(reg_filename.replace("_R_","_G_"), img_gs_reg_crop.astype(np.int16), plugin='tifffile', check_contrast=False)
        io.imsave(reg_filename.replace("_R_","_P_"), img_ps_reg_crop.astype(np.int16), plugin='tifffile', check_contrast=False)


# #### Segment the nuclear fluorescent images using cellpose  
# If mask files already exist, skip this step, otherwise:    
# Load the registered red stack nuclear marker images  
# Segment each image and only keep masks greater than a selected area  
# Save the mask files as an image sequence  
# For compatibility with the tracking method, save the masks and the nuclear intensity images as indivdual files  

# In[7]:


n_diameter = 13
flow_threshold = .7
mask_threshold=0
min_size=75
resample = True

#call cellpose on each image to segment the nuclei
# DEFINE CELLPOSE MODEL
model = models.Cellpose(gpu=True, model_type='cyto')

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# will use channel R = 1 as nuclear channel only
n_channels = [0,0]

for subdir in subdirectories:
    results = [] #collect results for one field in the well
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    mask_filename = reg_filename.replace("_reg_stack.tif", "_nuc_masks_stack.png")
    tracking_path = os.path.join(output_path,'tracking/')
    if not os.path.exists(tracking_path+well+"/"+field+"/nuc_masks/"):
            os.makedirs(tracking_path+well+"/"+field+"/nuc_masks/")
    if not os.path.exists(tracking_path+well+"/"+field+"/reg"):
            os.makedirs(tracking_path+well+"/"+field+"/reg")
    if not os.path.exists(tracking_path+well+"/"+field+"/results"):
            os.makedirs(tracking_path+well+"/"+field+"/results")
    if not os.path.exists(tracking_path+well+"/"+field+"/filtered_masks"):
            os.makedirs(tracking_path+well+"/"+field+"/filtered_masks")

    if not os.path.exists(mask_filename): #Only segment if no mask file
        print("Start segmenting "+reg_filename)
        img_rs_reg = io.imread(reg_filename)
        n_mask_images = []
        for i, image in enumerate(img_rs_reg):
            # create nuclear masks with cellpose 
            n_masks, flows, styles, diams = model.eval(image,
                                             diameter=n_diameter,
                                             flow_threshold=flow_threshold,
                                             mask_threshold=mask_threshold,
                                             channels=n_channels,
                                             min_size=min_size,
                                             resample = resample)
       
            #Save mask image
            n_mask_images.append(n_masks)
            
        tifffile.imwrite(mask_filename, np.array(n_mask_images), imagej=True)

        # Open the mask stack:
        im = Image.open(mask_filename)
 
        # create an index variable:
        i =0
        app = []
 
        # iterate over the mask stack and save each frame to disk:
        print("Save individual nuclear masks to disk for "+reg_filename)
        for fr in ImageSequence.Iterator(im):
            app.append(fr)
            fr.save(tracking_path+well+"/"+field+"/nuc_masks/"+"mask%03.d.tif"%i)
            i = i + 1
            
        # Open the registered intensity stack: use the phase images
        im = Image.open(reg_filename)
 
        # create an index variable:
        i = 0
        app = []
 
        # iterate over the registered phase stack and save each frame to disk:
        for fr in ImageSequence.Iterator(im):
            app.append(fr)
            fr.save(tracking_path+well+"/"+field+"/reg/"+"t%03.d.tif"%i)
            i = i + 1 


# #### Track the nuclei  
# Use the KIT-Loeffler tracking method to track the nuclei    
# The tracking output includes masks with new label values and a res_track.txt file as described below  

# In[ ]:


for subdir in subdirectories: #track all fileds in the well
    field = re.findall("field_[1-9]",subdir)[0]
    #Condition on the tracking output file res_tracks.txt existing
    if not os.path.exists(tracking_path+well+"/"+str(field)+"/results/res_track.txt"): #Only segment if no mask file
        cmd = "python -m run_tracking --image_path "+ tracking_path+well+"/"+str(field)+"/reg/ --segmentation_path "+tracking_path+well+"/"+str(field)+"/nuc_masks/ --results_path "+tracking_path+well+"/"+str(field)+"/results --delta_t 3 --default_roi_size 2"
        returned_value = os.system(cmd)  # returns the exit code in unix
        print('returned value from command line tracking with KIT:', returned_value)


# #### Identify cells
# Read in the tracking results  
# res_track.txt - A text file representing an acyclic graph for the whole image sequence. Every line corresponds to a single track that is encoded by four numbers separated by a space:  
# L B E P where  
# L - a unique label of the track (label of markers, 16-bit positive value)  
# B - a zero-based temporal index of the frame in which the track begins  
# E - a zero-based temporal index of the frame in which the track ends  
# P - label of the parent track (0 is used when no parent is defined)
# 
# 
# Filter the track objects keeping the parents and those with a minimum track length and save the results to tracks.csv    
# 
# Create a new file tracks.csv with the following columns:  
# label - a unique label of the track (label of markers, 16-bit positive value)  
# begins - a zero-based temporal index of the frame in which the track begins  
# ends - a zero-based temporal index of the frame in which the track ends  
# parent - label of the parent track (0 is used when no parent is defined)  
# length - The number of frames that the cell is identified in  
# plateID - Character string of the plate's ID such as AU02001  
# well - Character string of the well such as A1  
# field - Integer of the image field within the well  

# In[ ]:


#set filter parameters
min_track_length = 3
#loop through the results from each segmented field
for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    res_filename = os.path.join(output_path,"tracking",well,field,"results","res_track.txt")
    res_flt_filename = res_filename.replace("res_track.txt","tracks.csv")
    tracks = pd.read_csv(res_filename, sep=" ", header=None)
    tracks.columns = ["label", "begins", "ends", "parent"]
    tracks['length'] = tracks.ends - tracks.begins + 1
    last_track = tracks.ends.max()
    #check if object is a parent
    tracks["is_parent"] = tracks['label'].isin(tracks['parent'])
    tracks['plateID'] = plateID
    tracks['well'] = well
    tracks['field'] = field.replace("field_","")
    #If filtered results do not exist, read in the res_track.txt file for the current field
    if not os.path.exists(res_flt_filename):
        #Filter using the filter parameters
        #remove short tracks that are not parents and are not in the last frame
        tracks_flt = tracks.query('length >= @min_track_length or is_parent or ends == @last_track')
        ##remove any track that appears after the first frame and doesn't have a parent
        #tracks_flt = tracks_flt.query('not (begins > 1 & parent == 0)')
        #write out the res_flt_track.txt file
        tracks_flt.to_csv(res_flt_filename,index=False) 


# #### Filter masks to only tracked cells  
# Use the filtered tracks to remove masks for non-cell objects  
# Save the filtered masks as individual image files in filtered_masks directory   

# In[ ]:


#loop through the fields in the well
for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    res_flt_filename = os.path.join(output_path,"tracking",well,field,"results","tracks.csv")
    mask_track_path = os.path.join(output_path,"tracking",well,field,"results")
   #read in the tracks file for this field
    tracks = pd.read_csv(res_flt_filename) 
    #loop through the mask images in the field
    tracked_mask_filenames = sorted(glob.glob(mask_track_path+"/mask*"))
    # iterate over the mask files
    for fn in tracked_mask_filenames:
        #read in the mask image
        im = io.imread(fn)
        #replace any label that's not a cell with a 0 value
        cell_labels = np.array([x if x in tracks.label.to_numpy()
                                   else 0 for x in range(0, im.max()+1)])
        im_filtered = cell_labels[im]
        
        io.imsave(fn.replace("results","filtered_masks"), im_filtered.astype(np.int16), plugin='tifffile', check_contrast=False)
        


# #### Get excel metadata file  
# If this file does not exists, cretae a level 0 file that is data only  

# In[ ]:


#If the metadata exists, join it to the data and write out as a level 1 file
metadata_filename = os.path.join(data_path,plateID,"metadata",plateID+".xlsx")

if os.path.exists(metadata_filename):
    md_all = pd.read_excel(metadata_filename, engine='openpyxl', dtype={'Drug1Concentration': str, 'Drug2Concentration': str})
    
    #remove unwanted columns read in from the excel files
    r = re.compile("Unnamed.*")
    columns_to_drop = list(filter(r.match, md_all.columns)) 
    metadata = md_all.drop(columns = columns_to_drop)
    
    #match metadata and data well labels format
    metadata['row'] = [re.sub(r'[0-9]*', '', Well) for Well in metadata['Well']]
    metadata['column'] = [re.sub(r'[A-Z]', '', Well) for Well in metadata['Well']]
    metadata['column'] = [re.sub(r'\A0', '', row) for row in metadata['column']]
    metadata['well'] = metadata['row'] + metadata['column']
    


# #### Pull data from images  
# Apply the filtered masks to the registered red nuclear channel and recored each cells nuclear morphology, intensity and texture  
# Create cytoplasm masks by expanding the nuclear masks by a fixed amount or until they collide with another nuclear expansion  
# Use the nuclear and cytoplasmic masks to measure intensities in the green fluroescent images  
# Calculate intensity ratios between the cytomplasm and nuclei  
# If the metadata is available, merge it with the cell level data  
# Store the cell level data (and metadata) in a csv file where each row is a cell  
# Data feature values can be decoded as follows:  
# \<compartment>\_\<pipeline name>\_\<channel name>\_\<regionprops name>  
# compartment - Nuclei, Cyto or Cell  
# pipeline name - PC for python cellpose or other if added  
# channel name - NR for nuclear reporter, CC for cell cycle reporter or others if added  
# regionprops name - label passed through from the skimage measure.regionprops function https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops  
# 
# 

# In[ ]:


cyto_expansion = 5
minutes_between_images = 30

#loop through the fields in the well

for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    filtered_mask_path = os.path.join(output_path,"tracking",well,field,"filtered_masks")
    tracked_mask_filenames = sorted(glob.glob(filtered_mask_path+"/mask*"))
    l0_filename = os.path.join(data_path+plateID,"Analysis",pipeline_name,plateID+"_"+well+"_"+field+"_level_0.csv")
    img_gs_reg = io.imread(output_path+plateID+"_G_"+well+"_"+field.replace("field_","")+"_reg_stack.tif")
    # iterate over the mask files
    results = []
    for i, fn in enumerate(tracked_mask_filenames):
        #read in the mask image
        masks = io.imread(fn)
        #read in registered R images
        reg_fn = fn.replace("filtered_masks","reg")
        image = io.imread(reg_fn.replace("mask","t"))

        #measure reporter intensity and nuclear morphology, texture
        nuclei = measure.regionprops_table(masks, intensity_image=image,
                                           properties=('label',
                                                       'area','bbox_area','convex_area','centroid','eccentricity','equivalent_diameter','extent','feret_diameter_max','filled_area',
                                                        'major_axis_length','minor_axis_length','moments_hu','perimeter','perimeter_crofton','solidity',
                                                        'mean_intensity','max_intensity','min_intensity'))#expand the masks to get cytoplasmic regions
        nuclei_boundaries = segmentation.find_boundaries(masks, mode='thick')*masks
        nuclei_expansions = segmentation.expand_labels(masks, cyto_expansion) - masks + nuclei_boundaries
    
        # measure nuclear and cytoplasmic intensities and textures in the green channel
        nuclei_g = measure.regionprops_table(masks, intensity_image=img_gs_reg[i],
                                             properties=('label','centroid',
                                                         'mean_intensity','max_intensity','min_intensity'))
        nuclei_exp_g = measure.regionprops_table(nuclei_expansions, intensity_image=img_gs_reg[i],
                                                 properties=('label','centroid',
                                                             'mean_intensity','max_intensity','min_intensity'))
    
        # turn results into a dataframe
        nuclei_data = pd.DataFrame(nuclei)
        nuclei_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch1_name+'_'+col  for col in nuclei_data.columns if col not in ['label']}, inplace=True)
   
        nuclei_g_data = pd.DataFrame(nuclei_g)
        nuclei_g_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_g_data.columns if col not in ['label']}, inplace=True)
    

        nuclei_exp_g_data = pd.DataFrame(nuclei_exp_g)
        nuclei_exp_g_data.rename(columns={col: 'Cyto_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuclei_exp_g_data.columns if col not in ['label']}, inplace=True)
       
        # recover the well and field values and add them to the dataframe
        well = re.findall('/[A-Z][0-9]+/',reg_fn)[0]
        well = re.sub('/','', well)
        nuclei_data['well'] = well
        field = re.findall('field_[0-9]+',reg_fn)[0]
        field = int(re.sub('field_','', field))
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
        df = pd.merge(nuclei_data, nuclei_g_data, how="left", on=["label"])
        df_all = pd.merge(df, nuclei_exp_g_data, how="left", on=["label"])

        # append the dataframe to the results list
        results.append(df_all)
    
    #concatenate all of the results from the sequences in the well
    l0_image = pd.concat(results)

    #join with the tracking results to get lineage, parent frame lenght values
    tracks_filename = os.path.join(output_path,"tracking",well,"field_"+str(field),"results/tracks.csv")
    tracks = pd.read_csv(tracks_filename)
    l0 = pd.merge(l0_image, tracks, how="left", on=["label", "well", "field"])

    if os.path.exists(metadata_filename):
        #merge data and metadata on well values
        l1= pd.merge(l0, metadata, how="left", on=["well"]).round(decimals=2)
        l1.to_csv(l0_filename.replace('level_0','level_1'), index = False)
    else:
        print("no metadata file for "+plateID+" so creating level 0 file")
        l0 = l0.round(decimals=2)
        l0.to_csv(l0_filename, index = False)


# In[ ]:




