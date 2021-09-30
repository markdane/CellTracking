#!/usr/bin/env python
# coding: utf-8


#setup libraries
import pandas as pd
import numpy as np
import os, re, glob, shutil
from skimage import io, util
import numpy as np


#read in phase images and write back out as 16uint files
data_path = "../../../images/"

for plate_name in ["AU00702","AU00801", "AU00802","AU00901", "AU01001", "AU01002","AU01101", "AU01102"]:
    data_paths = glob.glob(os.path.join(data_path,plate_name,"*/*/*_P_*.tif"),recursive=True)
    print("reformatting phase files in "+plate_name, flush = True)
    for data_path in data_paths:
        image_8 = io.imread(data_path)
        image_16 = np.uint16(image_8)
        io.imsave(data_path, image_16, check_contrast=False)


