#!/bin/bash
                                                                                
echo "Run ilastik on multiple well_locations"
echo

for beacon in 42
do
  echo $i
  srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Pixel_Classifier.ilp --raw_data="Beacon-$beacon/Red_Reg/Scene1Interval*_TxRed.tif" --output_format="tif sequence" --output_filename_format=Beacon-$beacon/Red_Reg/Probabilities_{slice_index}.tif
  srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data="Beacon-$beacon/Red_Reg/Scene1Interval*_TxRed.tif" --prediction_maps="Beacon-$beacon/Red_Reg/Probabilities_*.tif" --export_source="Plugin" --export_plugin="CSV-Table"
  #srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data="Beacon-$beacon/Red_Reg/Scene1Interval*_TxRed.tif" --prediction_maps="Beacon-$beacon/Red_Reg/Probabilities_*.tif" --export_source="Tracking-Result"
  srun -c 8 /home/users/dane/ilastik-1.3.2b3-Linux/run_ilastik.sh --headless --project=HeLa_Tracking.ilp --raw_data="Beacon-$beacon/GFP_Reg/Scene1Interval*_GFP.tif" --prediction_maps="Beacon-$beacon/Red_Reg/Probabilities_*.tif" --export_source="Plugin" --export_plugin="CSV-Table"
done
