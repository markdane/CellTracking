#!/usr/bin/env python
# coding: utf-8

# ### Live Cell Imaging File Copy
# 
# This notebook intiates processing  image files by first copying them to the exacloud gscratch storage.

# In[2]:


#setup libraries
import numpy as np
import pandas as pd
import os, re, glob, shutil, sys


# In[ ]:


def get_file_df(plateID):
    print("getting files for "+plateID)
    #data_path =  '/home/groups/graylab_share/OMERO.rdsStore/liby/'+directory_name+'/'+plateID+"/"+plateID+'_IMAGES'
    #data_path =  '/home/groups/heiserlab_genomics/home/liby/'+directory_name+'/'+plateID+"/"+plateID+'_IMAGES'
    data_path =  '/home/groups/heiserlab_genomics/home/liby/'+directory_name+'/'+plateID+'_IMAGES' #Use for HCC1143 images
    #data_path = "/home/groups/heiserlab_genomics/home/grossse/AU565 Project/AU565 Image Files/AU_I_L_006_01_1/AU_I_L_006_01_1_IMAGES"
    #data_path =   '/home/groups/heiserlab_genomics/globus/calistri/'+directory_name+'/'+plateID

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
    
    print("get existing files")
    existing_file_paths = glob.glob(os.path.join(dest_path,plateID,"[A-D][1-6]","field_[1-4]","*.tif"),recursive=False)
    existing_df = pd.DataFrame(existing_file_paths, columns =['path'])
    existing_df['filename'] = existing_df['path'].str.extract('(\w*.tif)')
    existing_df['plateID'] = existing_df['filename'].str.extract('(^[a-zA-Z0-9]*)')
    existing_df['channel'] = existing_df['filename'].str.extract('(_[PRG]_)')
    existing_df['channel'] = existing_df['channel'].str.extract('([PRG])')
    existing_df['well'] = existing_df['filename'].str.extract('(_[A-Z][0-9]+_)')
    existing_df['well'] = existing_df['well'].str.extract('([A-Z][0-9]+)')
    existing_df['field'] = existing_df['filename'].str.extract('(_[0-9]+_)')
    existing_df['field'] = existing_df['field'].str.extract('([0-9]+)')
    existing_df['time'] = existing_df['filename'].str.extract('(_[a-z0-9]*.tif)')
    existing_df['time'] = existing_df['time'].str.replace('[_(.tif)]','',regex = True)
    #delete existing file paths from df to not waste time overwriting them
    #only execute if there are existing files
    if(existing_df.shape[0] >0):
        df = df.loc[np.isin(df['filename'].to_numpy(), existing_df['filename'].to_numpy(), assume_unique=True, invert= True)] 
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

def get_dest_df(dest_path):
    print("get existing files")
    existing_file_paths = glob.glob(os.path.join(dest_path,plate_name,"[A-D][1-6]","field_[1-4]","*.tif"),recursive=False)


dest_path = '/home/exacloud/gscratch/HeiserLab/images/'
directory_name = "AU565_drug_screen"
directory_name = "HCC1143_drug_screen"
#directory_name = "21MT1_drug_screen"
#directory_name = "MDAMB157_drug_screen"
#directory_name = "AU_I_L_006_01_1"
#directory_name = "nlc_incucyte"
#directory_name = "AU565_Drug_Screen_Project" #used for files stored on graylab_share

#plateID = "AU00601"
#plateID = "HC00701"
#plateID = "2101401"
#plateID = "MD00201"
#plateID = "221006_nlc_hcc1143ccrep_sirna_ptx_p2"

plateIDs = ["HC00701","HC00801","HC00901","HC01001","HC01301","HC01401"]
plateIDs = ["HC00901","HC01001","HC01301","HC01401"]
#plateIDs = ["2101001","2101201","2101301","2101401","2101501","2101601","2101701","2101801"]
#plateIDs = ["MD00301","MD00401","MD00501","MD00601","MD00701","MD00801"]

for plateID in plateIDs:
    print("starting file copy for "+plateID)
    df_all = get_file_df(plateID)
    #poor_quality_wells = ("A5", "A6", "B5", "B6","C5", "C6","D5", "D6")
    #poor_quality_wells = ("A1", "A2", "A3", "A4", "A5","A6",
    #                      "B1", "B2", "B3", "B4", "B5", "B6", 
    #                      "C1", "C2",       "C4", "C5", "C5",
    #                      "D1", "D2", "D3", "D4", "D5", "D6")

    poor_quality_wells = ()

    df = df_all[~df_all.well.isin(poor_quality_wells)].reset_index()
    df_directories = make_directories(df, dest_path)

    df = df.join(df_directories.set_index(['plateID','well','field']),on = ['plateID','well','field'])
    df['full_dest_path'] = dest_path+df['dest_path']+"/"+df['filename']
    df_copied = glob.glob(os.path.join(dest_path,plateID,"*/field*/*.tif"),recursive=True)
    df_copy = df[~df['full_dest_path'].isin(df_copied)]

    attempts = 0
    while (np.logical_and(len(df_copy) >0, attempts < 4)):
        print("copying files to "+dest_path+plateID)
        for count, src_path in enumerate(df_copy['path'].tolist()):
            shutil.copy2(src_path, df_copy.iloc[count]['full_dest_path'])
        #get files successfully copied
        df_copied = glob.glob(os.path.join(dest_path,plateID,"*/field*/*.tif"),recursive=True)
        #update dataframe of files that still need to be copied
        df_copy = df_copy[~df_copy['full_dest_path'].isin(df_copied)]
        attempts += 1
    print("done copying files for "+plateID+"\n")


# In[ ]:




