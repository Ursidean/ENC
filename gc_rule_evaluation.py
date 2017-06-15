"""
Evaluate a growing clusters rule set using Kappa, Kappa Simulation and clumpiness.
"""


import numpy as np
from considered_distances import considered_distances
import csv
from read_map import read_map
from kappa import kappa
from kappa import ksim
from clumpy_module import clumpiness_index
from run_metro import run_metro
from set_rand import set_rand
from set_NR import set_lp_rule


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\ENC\\"
project_file_basepath = "C:\\Users\\a1210607\\Geonamica\\Metronamica\\"
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
log_settings_2000 = base_path + "LogSettings.xml"
log_settings_2006 = base_path + "LogSettings2006.xml"
output_path = base_path + "Example_null_output\\"
base_seed = 1000
max_runs = 10

#Case study variables.
# Land_use_class_names
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures", "Agricultural areas",
             "Residential", "Industry & commerce", "Recreation areas", "Forest", "Road & rail",
             "Seaports", "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Maximum_neighbourhood_distance_considered
dmax = 5
# Number_of_land_use_classes
luc = len(luc_names)
# Number_of_passive_land_use_classes
pas = 1
# Number_of_feature_land_use_classes
fea = 6
# Number_of_active_land_use_classes
act = luc - (fea + pas)

# Data specification

case_studies = ["Berlin", "Budapest", "Madrid", "Rome"]
# Create empty dictionaries to store maps of data at 1990, 2000, and masks
maps_1990 = {}
maps_2000 = {}
maps_2006 = {}
masks = {}
data_clumpiness = {}
map_directory = (base_path + "Example_case_study_data\\")

# Read in data maps
for i in range(0, len(case_studies)):
    key = case_studies[i]
    map_1990_path = (map_directory + key + "\\" + key.lower() + "_1990.asc")
    map_2000_path = (map_directory + key + "\\" + key.lower() + "_2000.asc")
    map_2006_path = (map_directory + key + "\\" + key.lower() + "_2006.asc")
    mask_path = (map_directory + key + "\\" + key.lower() + "_mask.asc")
    maps_1990[key] = read_map(map_1990_path)
    maps_2000[key] = read_map(map_2000_path)
    maps_2006[key] = read_map(map_2006_path)
    masks[key] = read_map(mask_path)
    # Data results generation.

    # Specify omap, amap and mask
    map_1990 = read_map(map_1990_path)
    map_2000 = read_map(map_2000_path)
    map_2006 = read_map(map_2006_path)
    mask = read_map(mask_path)
    # Calculate data clumpiness,1990,2000,2006
    data_map_1990_clumpiness = clumpiness_index(map_1990, mask, luc)
    data_map_2000_clumpiness = clumpiness_index(map_2000, mask, luc)
    data_map_2006_clumpiness = clumpiness_index(map_2006, mask, luc)
    # Write to dictionary
    key1 = case_studies[i] + "-clumpiness-1990"
    data_clumpiness[key1] = data_map_1990_clumpiness
    key2 = case_studies[i] + "-clumpiness-2000"
    data_clumpiness[key2] = data_map_2000_clumpiness
    key3 = case_studies[i] + "-clumpiness-2006"
    data_clumpiness[key3] = data_map_2006_clumpiness

# Write output to a csv file

data_clumpiness_file = output_path + "data_clumpiness.csv"
store = [0] * (luc + 1)
with open(data_clumpiness_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Case study"] + luc_names
    writer.writerow(values)
    for i in range(0, len(case_studies)):
        key1 = case_studies[i] + "-clumpiness-1990"
        store[0] = key1
        for j in range(0, luc):
            store[j + 1] = data_clumpiness[key1][j]
        writer.writerow(store)
        key2 = case_studies[i] + "-clumpiness-2000"
        store[0] = key2
        for j in range(0, luc):
            store[j + 1] = data_clumpiness[key2][j]
        writer.writerow(store)
        key3 = case_studies[i] + "-clumpiness-2006"
        store[0] = key3
        for j in range(0, luc):
            store[j + 1] = data_clumpiness[key3][j]
        writer.writerow(store)

# Generate null model results
# Year 2000----------------------------------------------------------------------

# Create empty dictionary to store metric values
null_model_metrics = {}
for i in range(0, len(case_studies)):
    project_file_key = case_studies[i]

    # Specify case study variables
    project_file_path_2000 = project_file_basepath + case_studies[i]
    project_file_2000 = project_file_path_2000 + "\\" + case_studies[i] + ".geoproj"
    project_file_path_2006 = project_file_basepath + case_studies[i] + "_2006\\"
    project_file_2006 = project_file_path_2006 + case_studies[i] + "_2006.geoproj"

    # Generate null model rules, input
    rules = {}
    for p in range(0, luc):
        for q in range(0, act):
            if p == q + pas:
                key = "from " + luc_names[p] + " to " + luc_names[q + pas]
                rules[key] = [100, 1, 0, 5]
                fu_elem = q
                lu_elem = p
                y0 = rules[key][0]
                y1 = rules[key][1]
                y2 = rules[key][2]
                xe = rules[key][3]
                set_lp_rule(project_file_2000, fu_elem, lu_elem, y0, y1, y2, xe)
                set_lp_rule(project_file_2006, fu_elem, lu_elem, y0, y1, y2, xe)
            else:
                key = "from " + luc_names[p] + " to " + luc_names[q + pas]
                rules[key] = [0, 0, 0, 5]
                fu_elem = q
                lu_elem = p
                y0 = rules[key][0]
                y1 = rules[key][1]
                y2 = rules[key][2]
                xe = rules[key][3]
                set_lp_rule(project_file_2000, fu_elem, lu_elem, y0, y1, y2, xe)
                set_lp_rule(project_file_2006, fu_elem, lu_elem, y0, y1, y2, xe)
    # Specify the simulated map paths for 2000 & 2006
    smap_path_2000 = project_file_path_2000 + "\\Log\\Land_use\\Land use map_2000-Jan-01 00_00_00.rst"
    smap_path_2006 = project_file_path_2006 + "\\Log\\Land_use\\Land use map_2006-Jan-01 00_00_00.rst"
    # Specify land-use maps at 1990, 2000, 2006, and masking map
    map_1990 = maps_1990[project_file_key]
    map_2000 = maps_2000[project_file_key]
    map_2006 = maps_2006[project_file_key]
    mask = masks[project_file_key]
    # Generate simulated output maps, calculate results
    for c in range(0, max_runs):
        key = case_studies[i] + "-" + str(c)
        # Set random seed
        rseed = base_seed + c
        set_rand(project_file_2000, rseed)
        set_rand(project_file_2006, rseed)
        # Run metronamica to generate simulated output maps
        run_metro(project_file_2000, log_settings_2000, project_file_path_2000, geo_cmd)
        run_metro(project_file_2006, log_settings_2006, project_file_path_2006, geo_cmd)
        # Read simulated map
        null_model_2000_map = read_map(smap_path_2000)
        null_model_2006_map = read_map(smap_path_2006)
        # Calculate the metrics for the year 2000.
        it_kappa = kappa(map_2000, null_model_2000_map, mask)
        it_ksim = ksim(map_1990, map_2000, null_model_2000_map, mask)
        it_clu = clumpiness_index(null_model_2000_map, mask, luc)
        # Write to dictionary
        key1 = case_studies[i] + "-kappa-2000-it-" + str(c)
        null_model_metrics[key1] = it_kappa
        key2 = case_studies[i] + "-ksim-2000-it-" + str(c)
        null_model_metrics[key2] = it_ksim
        key3 = case_studies[i] + "-clu-2000-it-" + str(c)
        null_model_metrics[key3] = it_clu

        # Calculate the metrics for the year 2006.
        it_kappa = kappa(map_2006, null_model_2006_map, mask)
        it_ksim = ksim(map_2000, map_2006, null_model_2006_map, mask)
        it_clu = clumpiness_index(null_model_2006_map, mask, luc)
        # Write to dictionary
        key1 = case_studies[i] + "-kappa-2006-it-" + str(c)
        null_model_metrics[key1] = it_kappa
        key2 = case_studies[i] + "-ksim-2006-it-" + str(c)
        null_model_metrics[key2] = it_ksim
        key3 = case_studies[i] + "-clu-2006-it-" + str(c)
        null_model_metrics[key3] = it_clu

# Write output to csv file, 2000
null_models_2000_metrics_file = output_path + "null_model_metrics_2000.csv"
store = [0] * (4 + luc)
with open(null_models_2000_metrics_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Case study", "It.", "Kappa", "KSIM"] + luc_names
    writer.writerow(values)
    for i in range(0, len(case_studies)):
        for c in range(0, max_runs):
            key1 = case_studies[i] + "-kappa-2000-it-" + str(c)
            key2 = case_studies[i] + "-ksim-2000-it-" + str(c)
            key3 = case_studies[i] + "-clu-2000-it-" + str(c)
            store[0] = case_studies[i]
            store[1] = c
            store[2] = null_model_metrics[key1]
            store[3] = null_model_metrics[key2]
            for j in range(0, luc):
                store[4 + j] = null_model_metrics[key3][j]
            writer.writerow(store)

# Write output to csv file, 2006
null_models_2006_metrics_file = output_path + "null_model_metrics_2006.csv"
store = [0] * (4 + luc)
with open(null_models_2006_metrics_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Case study", "It.", "Kappa", "KSIM"] + luc_names
    writer.writerow(values)
    for i in range(0, len(case_studies)):
        for c in range(0, max_runs):
            key1 = case_studies[i] + "-kappa-2006-it-" + str(c)
            key2 = case_studies[i] + "-ksim-2006-it-" + str(c)
            key3 = case_studies[i] + "-clu-2006-it-" + str(c)
            store[0] = case_studies[i]
            store[1] = c
            store[2] = null_model_metrics[key1]
            store[3] = null_model_metrics[key2]
            for j in range(0, luc):
                store[4 + j] = null_model_metrics[key3][j]
            writer.writerow(store)