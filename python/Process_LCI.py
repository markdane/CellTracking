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
from scipy import stats, spatial
from PIL import Image, ImageSequence
import tifffile


# In[2]:


pipeline_name = "CK" #cellpose and KIT tracking
ch1_name = 'NR'
ch2_name = 'CC'
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
#data_path = '/Users/dane/Documents/CellTrackingProjects/AU565/images/'
plateID = 'AU02001'
plateID = sys.argv[1]
well_index = 1
well_index = int(sys.argv[2])

output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files/")
transformation_path = os.path.join(output_path,"transformations")
tracking_path = os.path.join(output_path,'tracking/')

#Only use subdirectories within the selected well
well_directory = sorted(glob.glob(os.path.join(data_path+plateID,"[A-Z][1-9]")))[well_index-1]
well = re.findall("[A-Z][1-9]$",well_directory)[0]
subdirectories = sorted(glob.glob(os.path.join(well_directory,"field_[1-9]")))[0:2] ######limit to two fields per well


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
 #       img_c_reg = np.stack((img_rs_reg.astype(np.int16), img_ps_reg.astype(np.int16)), axis = -1) 
#creates 21k files with (21, 1040, 1408, 2) shape
        img_c_reg = np.stack((img_rs_reg.astype(np.int16), img_ps_reg.astype(np.int16)), axis = 1) #create tcxy 


        if not os.path.exists(output_path):
            os.makedirs(output_path)

        print("saving "+reg_filename)
        io.imsave(reg_filename, img_rs_reg.astype(np.int16), plugin='tifffile', check_contrast=False)
        print("saving "+reg_filename.replace("_R_","_G_"))
        io.imsave(reg_filename.replace("_R_","_G_"), img_gs_reg.astype(np.int16), plugin='tifffile', check_contrast=False)
        print("saving "+reg_filename.replace("_R_","_P_"))
        io.imsave(reg_filename.replace("_R_","_P_"), img_ps_reg.astype(np.int16), plugin='tifffile', check_contrast=False)
        print("saving "+reg_filename.replace("_R_","_C_"))
        io.imsave(reg_filename.replace("_R_","_C_"), img_c_reg, plugin='tifffile', check_contrast=False)


# #### Segment the phase  images using cellpose  
# If mask files already exist, skip this step, otherwise:    
# Load the registered phase images  
# Segment each image and only keep masks greater than a selected area  
# Save the mask files as an image sequence  
# For compatibility with the tracking method, save the masks and the phase images as indivdual files  

# In[ ]:


diameter = 16
n_diameter = 13
flow_threshold = 0
n_flow_threshold = .7
mask_threshold=0
n_mask_threshold=0
min_size=75
n_min_size=75
resample = True

#call cellpose on each image to segment the phase cytoplasm and nuclei
# DEFINE CELLPOSE MODELs
model = models.Cellpose(gpu=True, model_type='cyto')
n_model = models.Cellpose(gpu=True, model_type='nuclei')

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# will use channel R = 1 as nuclear channel only
cy_channels = [2,1]
n_channels = [1,0]

for subdir in subdirectories:
    results = [] #collect results for one field in the well
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    reg_filename = os.path.join(output_path,plateID+"_C_"+well+"_"+field_num+"_reg_stack.tif")
    cy_mask_filename = reg_filename.replace("_reg_stack.tif", "_masks_stack.png")
    tracking_path = os.path.join(output_path,'tracking/')
    if not os.path.exists(tracking_path+well+"/"+field+"/masks/"):
            os.makedirs(tracking_path+well+"/"+field+"/masks/")
    if not os.path.exists(tracking_path+well+"/"+field+"/nmasks/"):
            os.makedirs(tracking_path+well+"/"+field+"/nmasks/")
    if not os.path.exists(tracking_path+well+"/"+field+"/reg"):
            os.makedirs(tracking_path+well+"/"+field+"/reg")
    if not os.path.exists(tracking_path+well+"/"+field+"/results"):
            os.makedirs(tracking_path+well+"/"+field+"/results")
    if not os.path.exists(tracking_path+well+"/"+field+"/filtered_masks"):
            os.makedirs(tracking_path+well+"/"+field+"/filtered_masks")

    if not os.path.exists(cy_mask_filename): #Only segment if no mask file
        img_rs_reg = io.imread(reg_filename)
                         
        cy_mask_images = []
        n_mask_images = []
        for i, image in enumerate(img_rs_reg):
            print("processing "+reg_filename+" index "+str(i))
            #Segment using the nuclear signal in the red channel
            #background = restoration.rolling_ball(image, radius=50)
            #image_bc = exposure.adjust_log(image-background, gain=1, inv=False)
            #image_rf = filters.median(image_bc, morphology.disk(1))

            # create masks with cellpose 
            cy_masks, flows, styles, diams = model.eval(image,
                                             diameter=diameter,
                                             flow_threshold=flow_threshold,
                                             mask_threshold=mask_threshold,
                                             channels=cy_channels,
                                             min_size=min_size,
                                             resample = resample)
       
            #Save mask image
            cy_mask_images.append(cy_masks)
            
             # create nuclear masks with cellpose 
            n_masks, flows, styles, diams = n_model.eval(image,
                                             diameter=n_diameter,
                                             flow_threshold=n_flow_threshold,
                                             mask_threshold=n_mask_threshold,
                                             channels=n_channels,
                                             min_size=n_min_size,
                                             resample = resample)
       
            #Save mask image
            n_mask_images.append(n_masks)
            
        tifffile.imwrite(cy_mask_filename, np.array(cy_mask_images), imagej=True)
        tifffile.imwrite(cy_mask_filename.replace('masks', 'nmasks'), np.array(n_mask_images), imagej=True)

        # Open the cyto mask stack:
        im = Image.open(cy_mask_filename)
 
        # create an index variable:
        i =0
        app = []
 
        # iterate over the cyto mask stack and save each frame to disk:
        for fr in ImageSequence.Iterator(im):
            app.append(fr)
            fr.save(tracking_path+well+"/"+field+"/masks/"+"mask%03.d.tif"%i)
            i = i + 1
            
         # Open the nuclear mask stack:
        im = Image.open(cy_mask_filename.replace('masks', 'nmasks'))
 
        # create an index variable:
        i =0
        app = []
 
        # iterate over the nuclear mask stack and save each frame to disk:
        for fr in ImageSequence.Iterator(im):
            app.append(fr)
            fr.save(tracking_path+well+"/"+field+"/nmasks/"+"nmask%03.d.tif"%i)
            i = i + 1
        
#        # Open the registered intensity stack: use the phase images
        im = Image.open(reg_filename.replace("_C_", "_P_"))
 
        # create an index variable:
        i = 0
        app = []
 
        # iterate over the registered phase stack and save each frame to disk:
        for fr in ImageSequence.Iterator(im):
            app.append(fr)
            fr.save(tracking_path+well+"/"+field+"/reg/"+"t%03.d.tif"%i)
            i = i + 1 


# #### Track the nuclei  
# Use the KIT-Loeffler tracking method to track the masks and the nuclei    
# For now, use the string below to run from the tracking script from the command line  
# The tracking output includes masks with new label values and a res_track.txt file as described below  

# In[ ]:


cmd = "python -m run_tracking --image_path "+ tracking_path+well+"/"+str(field)+"/reg/ --segmentation_path "+tracking_path+well+"/"+str(field)+"/masks/ --results_path "+tracking_path+well+"/"+str(field)+"/results --delta_t 3 --default_roi_size 2"
returned_value = os.system(cmd)  # returns the exit code in unix
print('returned value:', returned_value)


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
    #check if object is a parent
    tracks["is_parent"] = tracks['label'].isin(tracks['parent'])
    tracks['plateID'] = plateID
    tracks['well'] = well
    tracks['field'] = field.replace("field_","")
    #If filtered results do not exist, read in the res_track.txt file for the current field
    if not os.path.exists(res_flt_filename):
        #Filter using the filter parameters
        #remove short tracks that are not parents
        tracks_flt = tracks.query('length >= @min_track_length or is_parent')
        ##remove any track that appears after the first frame and doesn't have a parent
        #tracks_flt = tracks_flt.query('not (begins > 1 & parent == 0)')
        #write out the res_flt_track.txt file
        tracks_flt.to_csv(res_flt_filename,index=False) 


# #### Filter masks to only tracked cells  
# Use the filtered tracks to remove masks for non-cell objects  
# Save the filtered masks as individual image files in filtered_masks directory   

# In[ ]:


#TODO make conditional on no filtered mask files#loop through the fields in the well
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
        


# #### Find nuclei within each cytoplasmic mask  
# Index through each filtered cytoplasmic mask and apply it to the red nuclear channel. Then run cellpose on the results to get nuclei masks. Count the number of masks in the cell and select one to represent the cell's nuclei. Either rebuild a single nuclear mask image to be applied alongside the cytoplasmic mask or pull the region properties and build a dataframe with enough metadata to join it with the cytoplasmic region data. Save the nuclear mask image and data.  
# 

# ndiameter = 6
# flow_threshold = 2
# mask_threshold=0
# min_size=10
# resample = True
# 
# #call cellpose on each image to segment the nuclei
# # DEFINE CELLPOSE MODEL
# model = models.Cellpose(gpu=True, model_type='nuclei')
# 
# # define CHANNELS to run segementation on
# # grayscale=0, R=1, G=2, B=3
# # channels = [cytoplasm, nucleus]
# # if NUCLEUS channel does not exist, set the second channel to 0
# # will use grayscale image as nuclear channel only
# nchannels = [0,0]
# 
# #loop through all fields in the well
# for subdir in subdirectories:
#     field = re.findall("field_[1-9]",subdir)[0]
#     nmask_filtered_path = os.path.join(output_path,"tracking",well,field,"nmasks")
#     nmask_filtered_filenames = glob.glob(nmask_filtered_path+"/nmask*") #see if any files exist
# if len(nmask_filtered_filenames) == 0: #Only segment if no nuclear mask file
#     print("creating nuclear masks for "+plateID+" "+well+" "+field)
#     #loop through the cyto mask images in the field
#     cmask_filenames = sorted(glob.glob(nmask_filtered_path.replace("nmasks","filtered_masks")+"/mask*"))
#     img_rs_reg = io.imread(output_path+plateID+"_R_"+well+"_"+field.replace("field_","")+"_reg_stack.tif")
# 
#     for i, fn in enumerate(cmask_filenames[0:2]):
#         print("reading cyto masks " + fn)
#         im_cmask = io.imread(fn)
#         #step through each non-zero label
#         for cy_label in np.unique(im_cmask)[200:201]:
#             if not cy_label == 0:
#                 print("processing label " + str(cy_label))
#                 label_mask = (im_cmask == cy_label) #create a full size mask with true values only for the current label
#                 label_mask = np.ones(img_rs_reg[i].shape)*label_mask
#                 nuc_int = np.where(im_cmask == cy_label,img_rs_reg[i], np.zeros(img_rs_reg[i].shape))
#                 #nuc_int = img_rs_reg[i]*label_mask #apply the full size mask and get a full size intensity image with only intensties within the mask area
#                 #run cellpose nuclei model on the image
#                 nmasks, flows, styles, diams = model.eval(nuc_int,
#                                              diameter=ndiameter,
#                                              flow_threshold=flow_threshold,
#                                              mask_threshold=mask_threshold,
#                                              channels=nchannels,
#                                              min_size=10,
#                                              resample = resample)
#                 
# 

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
# Apply the filtered masks to the registered cytoplasmic channel and record each cell's cytoplasmic morphology, intensity and texture  
# Apply the nuclear masks to the registered red channel and record each cell's nuclear morphology, intensity and texture  
# If the metadata is available, merge it with the cyoplasmic data   
# Store the cytoplasmic cell level data (and metadata) in a csv file where each row is a cell  
# Store the nuclear cell level data (and limited metadata) in a csv file where each row is a cell  
# Data feature values can be decoded as follows:  
# \<compartment>\_\<pipeline name>\_\<channel name>\_\<regionprops name>  
# compartment - Nuclei, Cyto or Cell  
# pipeline name - CK for cellpose KIT or other if added  
# channel name - NR for nuclear reporter, CC for cell cycle reporter or others if added  
# regionprops name - label passed through from the skimage measure.regionprops function https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops  
# 
# 

# In[ ]:


minutes_between_images = 30

#loop through the fields in the well

for subdir in subdirectories:
    field = re.findall("field_[1-9]",subdir)[0]
    filtered_mask_path = os.path.join(output_path,"tracking",well,field,"filtered_masks")
    cy_tracked_mask_filenames = sorted(glob.glob(filtered_mask_path+"/mask*"))
    l0_filename = os.path.join(data_path+plateID,"Analysis",pipeline_name,plateID+"_"+well+"_"+field+"_cell_level_0.csv")
    img_ps_reg = io.imread(output_path+plateID+"_P_"+well+"_"+field.replace("field_","")+"_reg_stack.tif")
    img_rs_reg = io.imread(output_path+plateID+"_R_"+well+"_"+field.replace("field_","")+"_reg_stack.tif")
    img_gs_reg = io.imread(output_path+plateID+"_G_"+well+"_"+field.replace("field_","")+"_reg_stack.tif")
    # iterate over the cyto mask files
    cell_results = []
    for i, fn in enumerate(cy_tracked_mask_filenames):
        #read in the cytoplasm mask image
        masks = io.imread(fn)
 
        #measure cell cycle reporter intensity and cellular morphology, texture
        cyto_cell_cycle = measure.regionprops_table(masks, intensity_image=img_gs_reg[i],
                                           properties=('label',
                                                       'area','bbox_area','convex_area','centroid','eccentricity','equivalent_diameter','extent','feret_diameter_max','filled_area',
                                                        'major_axis_length','minor_axis_length','moments_hu','perimeter','perimeter_crofton','solidity',
                                                        'mean_intensity','max_intensity','min_intensity'))
    
        # turn results into a dataframe
        cyto_data = pd.DataFrame(cyto_cell_cycle)
        cyto_data.rename(columns={col: 'Cell_'+pipeline_name+'_' +ch2_name+'_'+col  for col in cyto_data.columns if col not in ['label']}, inplace=True)
       
        # recover the well, field and time slice values and add them to the dataframe
        well = re.findall('/[A-Z][0-9]+/',fn)[0]
        well = re.sub('/','', well)
        cyto_data['well'] = well
        field = re.findall('field_[0-9]+',fn)[0]
        field = int(re.sub('field_','', field))
        cyto_data['field'] = field
        cyto_data['slice'] = i
        elapsed_minutes = i*minutes_between_images #assumes time slice numbering starts at 1
        day = np.floor(elapsed_minutes/(24*60)).astype(int)
        hour = np.floor((elapsed_minutes-day*(24*60))/60).astype(int)
        minute = np.floor(elapsed_minutes-day*(24*60)-hour*60).astype(int)
        day = str(day).zfill(2)
        hour = str(hour).zfill(2)
        minute = str(minute).zfill(2)
        cyto_data['time_slice'] = day+"d"+hour+"h"+minute+"m"
        
        #read in the nuclear mask image
        nfn = fn.replace("filtered_masks","nmasks")
        nfn = nfn.replace("/mask","/nmask")
        nmasks = io.imread(nfn)
 
        #measure cell cycle reporter intensity and nuclear morphology, texture
        nuc_cell_cycle = measure.regionprops_table(nmasks, intensity_image=img_gs_reg[i],
                                           properties=('label',
                                                       'area','bbox_area','convex_area','centroid','eccentricity','equivalent_diameter','extent','feret_diameter_max','filled_area',
                                                        'major_axis_length','minor_axis_length','moments_hu','perimeter','perimeter_crofton','solidity',
                                                        'mean_intensity','max_intensity','min_intensity'))
    
        # turn results into a dataframe
        nuc_data = pd.DataFrame(nuc_cell_cycle)
        nuc_data.rename(columns={col: 'Nuclei_'+pipeline_name+'_' +ch2_name+'_'+col  for col in nuc_data.columns if col not in ['label']}, inplace=True)
        nuc_data.rename(columns={"label": "nuc_label"}, inplace=True)
        # add the well, field and time slice values and add them to the dataframe
        #nuc_data['well'] = cyto_data['well']
        #nuc_data['field'] = cyto_data['field']
        #nuc_data['slice'] = cyto_data['slice']
        #nuc_data['time_slice'] = cyto_data['time_slice']

        #compute the distances between the cyto and nuclei centroids
        #Rows are cyto and columns are nuclei
        #smallest value in each row is the nearest nuclei to the cyto mask's centroid
        centroid_distances = pd.DataFrame(spatial.distance_matrix(cyto_data[['Cell_CK_CC_centroid-1', 'Cell_CK_CC_centroid-0']].to_numpy(),
        nuc_data[['Nuclei_CK_CC_centroid-1', 'Nuclei_CK_CC_centroid-0']].to_numpy()))

        #for each row, find the column with the smallest value
        nearest_nuclei = centroid_distances.idxmin(axis = 1)
        
        #add the nuclear label to the cyto data
        cyto_data['nuclei_label'] = nearest_nuclei
        
        #merge the cyto and nuclear data on the label
        cell_data = pd.merge(cyto_data, nuc_data, left_on='nuclei_label', right_on='nuc_label')
        
        # append the dataframe to the results list
        cell_results.append(cell_data)

    #concatenate all of the results from the images in the field
    l0 = pd.concat(cell_results)

    #join with the tracking results to get lineage, parent, frame length values
    tracks_filename = os.path.join(output_path,"tracking",well,"field_"+str(field),"results/tracks.csv")
    tracks = pd.read_csv(tracks_filename)
    l0 = pd.merge(l0, tracks, how="left", on=["label", "well", "field"])

    if os.path.exists(metadata_filename):
        #merge data and metadata on well values
        l1= pd.merge(l0, metadata, how="left", on=["well"]).round(decimals=2)
        l1.to_csv(l0_filename.replace('level_0','level_1'), index = False)
    else:
        print("no metadata file for "+plateID+" so creating level 0 file")
        l0 = l0.round(decimals=2)
        l0.to_csv(l0_filename, index = False)
    


# In[ ]:





# In[ ]:




