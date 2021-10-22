#!/usr/bin/env python
# coding: utf-8


import os, re, glob, sys, time

pipeline_name = "PI"
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
plateID = sys.argv[1]
output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
file_paths = sorted(glob.glob(output_path+"/*stack.tif"))

for file_path in file_paths:
    mask_file_path = re.sub('.tif','_Simple Segmentation.h5',file_path)
    well = re.findall("_[A-Z][1-9]_",mask_file_path)[0]
    well = re.sub("_","",well)
    field = re.findall("_[1-9]_",mask_file_path)[0]
    field = re.sub("_","", field)[0]
    if not os.path.exists(mask_file_path):
        cmd = "srun -c 4 --job-name=i"+plateID[3:8]+well+field+" -o i"+plateID[3:8]+well+field+".txt /home/users/dane/ilastik-1.4.0b15-Linux/run_ilastik.sh --headless --readonly=True --input_axes=tyxc --project=/home/groups/heiserlab_genomics/home/dane/CellTracking/AU565/AU00901_Pixel_112920_TRAINING.ilp --raw_data="+file_path+" --export_source='Simple Segmentation' --output_format='hdf5' --output_filename_format="+output_path+"/{nickname}_{result_type} &"
        returned_value = os.system(cmd)  # returns the exit code in unix
        if returned_value == 0:
            print("launched job to create masks for "+file_path)
        else:
            print("failed to launch job to create masks for "+file_path+" returned value "+str(returned_value))
        time.sleep(.1)

