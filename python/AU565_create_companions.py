#!/usr/bin/env python
# coding: utf-8

from ome_model.experimental import Plate, Image, create_companion
import subprocess
import re, sys
import pandas as pd

data_path = "AU565_library.csv"
plate_name = sys.argv[1]
file_metadata = pd.read_csv(data_path, dtype={'column': 'uint8'})
file_metadata = file_metadata[file_metadata['plate_name']==plate_name]

columns = pd.unique(file_metadata['column'])
rows = pd.unique(file_metadata['row'])
fields = pd.unique(file_metadata['field'])
timepoints = pd.unique(file_metadata['timepoint'])
channels = pd.unique(file_metadata['channel_name'])

channels = ['G', 'R', 'P']
PATTERN = re.compile(r"(\d{2})d(\d{2})h(\d{2})m")
print(f"Creating {plate_name}.companion.ome ...")


plate = Plate(plate_name, len(rows), len(columns))
for row_index, row in enumerate(rows):
    for column_index, column in enumerate(columns):
        well = plate.add_well(row_index, column_index)
        for field_index, field in enumerate(fields):
            # Create multi-channel timelapse image per field of view
            image_name = f"{plate_name}_{row}{column}_{field}"
            image = Image(
                image_name, 1408, 1040, 1, len(channels), len(timepoints),
                order="XYZTC", type="uint16")
            # Create channels
            image.add_channel(
                name='G',
                color = int.from_bytes([0, 255, 0, 255], 'big', signed=True),
                samplesPerPixel=1)
            image.add_channel(
                name='R',
                color = int.from_bytes([255, 0, 0, 255], 'big', signed=True),
                samplesPerPixel=1)
            image.add_channel(
                name='P',
                color = int.from_bytes([255, 255, 255, 255], 'big', signed=True),
                samplesPerPixel=1)
            tiff_folder = f"{plate_name}/{row}{column}/field_{field}/"
            for t, timepoint in enumerate(timepoints):
                # Parse timestamp from filename
                m = PATTERN.match(timepoint)
                t_min = (
                    int(m.group(1)) * 60 * 24 +
                    int(m.group(2)) * 60 +
                    int(m.group(3)))
                for c, channel in enumerate(channels):
                    # Define TIFF files for each plane
                    tiff_filename = (f"{plate_name}_{channel}_{row}{column}"
                                     f"_{field}_{timepoint}.tif")
                    image.add_tiff(tiff_folder + tiff_filename, c=c, z=0, t=t)
                    # Populate pixel size (in microns)
                    image.data['Pixels']['PhysicalSizeX'] = '1.24'
                    image.data['Pixels']['PhysicalSizeY'] = '1.24'
                    # Add plane metadata for each plane with timestamp
                    options = {
                        'DeltaT': f"{t_min}",
                        'DeltaTUnit': 'min',
                    }
                    image.add_plane(c=c, z=0, t=t, options=options)
            well.add_wellsample(field_index, image)


companion_file = f"{plate_name}.companion.ome"
create_companion(plates=[plate], out=companion_file)

# Indent XML for readability
proc = subprocess.Popen(
    ['xmllint', '--format', '-o', companion_file, companion_file],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE)
(output, error_output) = proc.communicate()

print("Done.")


