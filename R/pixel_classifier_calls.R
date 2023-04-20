#Call ilastik from the command line 

# lapply(6, function(beacon){
#   #paste0('/Applications/ilastik-1.3.2b3-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=HeLa_Pixels_Classifier.ilp --output_format="hdf5" --output_filename_format=\"Beacon-',beacon,'/Red_Reg/results/{nickname}_{slice_index}_results.h5\" \"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\"')
#   system(paste0('/Applications/ilastik-1.3.2b3-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=HeLa_Pixels_Classifier.ilp --output_format=hdf5 --output_filename_format=\"Beacon-',beacon,'/Red_Reg/results/Beacon_',beacon,'_results.h5\" \"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\"'))
# })

lapply(1, function(beacon){
  system(paste0('/Applications/ilastik-1.3.2b3-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=HeLa_Pixel_Classifier.ilp --raw_data=\"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\" --output_format=\"tif sequence\" --output_filename_format=Beacon-',beacon,'/Red_Reg/Probabilities_{slice_index}.tif'))
  system(paste0('/Applications/ilastik-1.3.2b3-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data=\"Beacon-',beacon,'/Red_Reg/Scene1Interval*_TxRed.tif\" --prediction_maps=\"Beacon-',beacon,'/Red_Reg/Probabilities_*.tif\" --export_source="Plugin" --export_plugin="CSV-Table"'))
  system(paste0('/Applications/ilastik-1.3.2b3-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data=\"Beacon-',beacon,'/GFP_Reg/Scene1Interval*_GFP.tif\" --prediction_maps=\"Beacon-',beacon,'/Red_Reg/Probabilities_*.tif\" --export_source="Plugin" --export_plugin="CSV-Table"'))
})
