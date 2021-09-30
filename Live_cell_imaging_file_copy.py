#!/usr/bin/env python
# coding: utf-8

# ### Live Cell Imaging File Copy
# 
# This script intiates processing AU565 image files by first copying them to the exacloud gscratch storage.


#setup libraries
import pandas as pd
import numpy as np
import os, re, glob, shutil, sys


plate_name = sys.argv[2]
directory_name = sys.argv[1]
#directory_name = "AU_I_L_015_02_1"
#plate_name = "AU01502"

def get_file_df(plate_name):
    print("getting files for "+plate_name)
    data_path =  '/home/groups/heiserlab_genomics/home/grossse/'+directory_name+'/'+plate_name+'_IMAGES'
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

print("starting file copy for "+plate_name)
df = get_file_df(plate_name)
df_directories = make_directories(df, dest_path)
    
df = df.join(df_directories.set_index(['plateID','well','field']),on = ['plateID','well','field'])
print("copying files to "+dest_path+'/'+plate_name)
for count, src_path in enumerate(df['path'].tolist()):
    full_dest_path = os.path.join(dest_path, df['dest_path'][count])
    shutil.copy2(src_path, full_dest_path)
print("done copying files for "+plate_name+"\n")