#!/usr/bin/env python
# coding: utf-8

# In[35]:


import os, re, glob, sys


# In[36]:


pipeline_name = 'PI' # python + ilastik
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
plateID = sys.argv[1]
input_files_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
image_stack_paths = sorted(glob.glob(input_files_path+"/*stack.tif"))
mask_paths = sorted(glob.glob(input_files_path+"/*Segmentation.h5"))
mainpath = os.path.join(data_path,plateID,"Analysis",pipeline_name)


# Only create an  job if there is no output file.

# In[59]:


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


# In[ ]:




