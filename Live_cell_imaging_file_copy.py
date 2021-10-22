#!/usr/bin/env python
# coding: utf-8

# ### Live Cell Imaging File Copy
# 
# This script intiates processing AU565 image files by first copying them to the exacloud gscratch storage.


#setup libraries
import pandas as pd
import numpy as np
import os, re, glob, shutil, sys
from skimage import io


plateID = sys.argv[2]
directory_name = sys.argv[1]
#directory_name = "AU_I_L_015_02_1"
#plateID = "AU01502"

def get_file_df(plateID):
    print("getting files for "+plateID)
    data_path =  '/home/groups/graylab_share/OMERO.rdsStore/liby/AU565_Drug_Screen_Project/'+directory_name+'/'+plateID+'_IMAGES'
    data_paths = glob.glob(os.path.join(data_path,"*.tif"),recursive=False)
    df = pd.DataFrame(data_paths, columns =['path'])
    df['filename'] = df['path'].str.extract('(\w*.tif)')
    df['plateID'] = df['filename'].str.extract('(^[a-zA-Z0-9]*)')
    df['channel'] = df['filename'].str.extract('(_[PRG]_)')
    df['channel'] = df['channel'].str.extract('([PRG])')
    df['well'] = df['filename'].str.extract('(_[A-Z][0-9]+_)')
    df['well'] = df['well'].str.extract('([A-Z][0-9]+)')
    df['field'] = df['filename'].str.extract('(_[0-9]+_)')
    df['field'] = df['field'].str.extract('([0-9]+)')
    df['time'] = df['filename'].str.extract('(_[a-z0-9]*.tif)')
    df['time'] = df['time'].str.replace('[_(.tif)]','',regex = True)
    return df

def make_directories(df, dest_path):
    print("making directories in "+dest_path)
    df_directories = df[['plateID', 'well', 'field']].drop_duplicates()
    df_directories['dest_path'] = df_directories['plateID']+'/'+df_directories['well']+'/'+'field_'+df_directories['field']
    
    for dest in df_directories['dest_path'].tolist():
        full_path = os.path.join(dest_path, dest)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
    return df_directories 
        
dest_path = '/home/exacloud/gscratch/HeiserLab/images/'

print("starting file copy for "+plateID)
df = get_file_df(plateID)
df_directories = make_directories(df, dest_path)
    
df = df.join(df_directories.set_index(['plateID','well','field']),on = ['plateID','well','field'])
print("copying files to "+dest_path+'/'+plateID)
for count, src_path in enumerate(df['path'].tolist()):
    full_dest_path = os.path.join(dest_path, df['dest_path'][count])
    shutil.copy2(src_path, full_dest_path)
print("done copying files for "+plateID+"\n")

pipeline_name = "PI" #python and ilastik
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
subdirectories = sorted(glob.glob(os.path.join(data_path+plateID,"[A-Z][1-9]","field_[1-9]")))

# Only process the files if there is no stack of RGP files

flourescent_scaler = 255/4095 #rescale from 12 to 8 bits
for subdir in subdirectories:
    well = re.findall("/[A-Z][1-9]",subdir)[0]
    well = re.findall("[A-Z][1-9]", well)[0]
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    c_filename = os.path.join(output_path,plateID+"_RGP_"+well+"_"+field_num+"_stack.tif")
    metadata_filename = os.path.join(output_path,plateID+"_RGP_"+well+"_"+field_num+"_filenames.txt")
    if not os.path.exists(c_filename):
        print("combining channels in "+subdir)
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

        img_c = np.stack([img_rs.astype('B'), img_gs.astype('B'), img_ps.astype('B')], axis = -1) 
    
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        io.imsave(c_filename, img_c, plugin='tifffile')
        with open(metadata_filename, 'w') as filehandle:
            for listitem in img_r_ic.files:
                filehandle.write('%s\n' % re.sub("_R_","_RGP_",listitem))

input_files_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
image_stack_paths = sorted(glob.glob(input_files_path+"/*stack.tif"))
mask_paths = sorted(glob.glob(input_files_path+"/*Segmentation.h5"))
mainpath = os.path.join(data_path,plateID,"Analysis",pipeline_name)

# Only create an  job if there is no output file.

wells = sorted(set(re.findall(r"_[A-Z][1-9]_", ''.join(mask_paths)))) #get a unique set of the wells with images that have been segmented

for well in wells:
    #initate a job if the csv file for the current well does not exists
    l1_file_path = data_path+plateID+"/Analysis/"+pipeline_name+"/"+plateID+well+"level_1.csv"
    well = re.sub("_","",well)
    if not os.path.exists(l1_file_path):
        cmd = 'srun -c 8 -J M'+plateID[3:7]+well[0:2]+' -o M'+plateID[3:7]+well[0:2]+'_out.txt -t 23:00:00 python Apply_ilastik_masks.py '+plateID+' '+well+' &'
        print("launching job with command "+cmd)
        returned_value = os.system(cmd)  # returns the exit code in unix
        if returned_value == 0:
            print("launched job to create data for "+l1_file_path)
        else:
            print("failed to launch job to create data for "+l1_file_path)


