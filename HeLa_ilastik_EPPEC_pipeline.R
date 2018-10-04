#ssh eppec.ohsu.edu
#cd /graylab/share/dane/HeLa_1
#srun -c 8 Rscript HeLa_ilastik_EPPEC_pipeline.R
#####DEBUG not working yet... may be due to r library upgrade


lapply(2, function(beacon){
  system(paste0('/home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Pixel_Classifier.ilp --raw_data=\"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\" --output_format=\"tif sequence\" --output_filename_format=Beacon-',beacon,'/Red_Reg/Probabilities_{slice_index}.tif'))
  system(paste0('/home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data=\"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\" --prediction_maps=\"Beacon-',beacon,'/Red_Reg/Probabilities_*.tif\" --export_source="Plugin" --export_plugin="CSV-Table"'))
  system(paste0('/home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data=\"Beacon-',beacon,'/GFP_Reg/Scene1Interval*_GFP.tif\" --prediction_maps=\"Beacon-',beacon,'/Red_Reg/Probabilities_*.tif\" --export_source="Plugin" --export_plugin="CSV-Table"'))
})
