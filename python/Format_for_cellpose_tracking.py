#!/usr/bin/env python
# coding: utf-8

#setup libraries
import numpy as np
import os, re, glob, sys
from skimage import io
from pystackreg import StackReg

pipeline_name = "PC" #python and cellpose
data_path = '/home/exacloud/gscratch/HeiserLab/images/'
plateID = sys.argv[1]

output_path = os.path.join(data_path+plateID,"Analysis",pipeline_name,"intermediate_files")
subdirectories = sorted(glob.glob(os.path.join(data_path+plateID,"[A-Z][1-9]","field_[1-9]")))

flourescent_scaler = 255/4095 #rescale from 12 to 8 bits

for subdir in subdirectories:
    well = re.findall("/[A-Z][1-9]",subdir)[0]
    well = re.findall("[A-Z][1-9]", well)[0]
    field = re.findall("field_[1-9]",subdir)[0]
    field_num = re.findall("[0-9]", field)[0]
    reg_filename = os.path.join(output_path,plateID+"_R_"+well+"_"+field_num+"_reg_stack.tif")
    # Only process the field-level image files if there is no stack of RGP files
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
        np.save(plateID+"_"+well+"_"+field+"_transformation_matrices.npy", tmats)
        
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

