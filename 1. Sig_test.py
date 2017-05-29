# This program conducts the significance test on the specified data.
import gdal
import numpy as np
import math
import csv
import time

# Specify the base path to the directory containing the empirical neighbourhood calibration tool-pack.
base_path = "C:\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Berlin"
# Set the paths to the directories and relevant data
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures", "Agricultural areas",
             "Residential", "Industry & commerce", "Recreation areas", "Forest", "Road & rail",
             "Seaports", "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive, feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5

# Read in the map for the data at time slice 1.
src_map = gdal.Open(map1_path)
omap = np.array(src_map.GetRasterBand(1).ReadAsArray())
# Read in the map for the data at time slice 2.
src_map = gdal.Open(map2_path)
amap = np.array(src_map.GetRasterBand(1).ReadAsArray())
# Read in the masking map.
src_map = gdal.Open(mask_path)
mask = np.array(src_map.GetRasterBand(1).ReadAsArray())

# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Determine the distances that will be analysed, use module: considered_distances.
from considered_distances import considered_distances
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
