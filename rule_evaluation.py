"""
Evaluate the input rule set using Kappa, Kappa Simulation and clumpiness.
"""


import numpy as np
from considered_distances import considered_distances
import csv        
from read_map import read_map
from kappa import kappa
from kappa import ksim
from clumpy_module import clumpiness_index
from area_weighted_clu import area_weighted_clu_error
from run_metro import run_metro
from set_rand import set_rand
from set_NR import set_lp_rule


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Berlin"
# Set the paths to the relevant directories
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
# Set the paths to the relevant files & executables
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
map3_path = data_path + case_study + "\\" + case_study.lower() + "_2006.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Seaports",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5
# Specify the fine-tuning calibration method parameters, the base random seed, the
# maximum number of simulation runs, and the golden section search tolerance.
max_runs = 50
base_seed = 1000
print "Testing case study:" + case_study
# Read in the maps and mask.
map1 = read_map(map1_path)
map2 = read_map(map2_path)
map3 = read_map(map3_path)
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(map1)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Determine the distances that will be analysed, use module: considered_distances.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])
# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[map2[i, j]] = luc_count[map2[i, j]] + 1

# Input the rules to be analysed. Start by initialising a dictionary for
# storage.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, 5]
# Load the attraction rules file
att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
att_rules = np.loadtxt(att_rule_file)
# Load the rules file, reading the inputs from the specified .csv file.
rules_file = output_path + case_study + "\\Rules\\final_rules.csv"
with open(rules_file, 'rb') as f:
    readCSV = csv.reader(f)
    next(f)
    for row in readCSV:
        i = row[0]
        j = row[1]
        key = "from " + i + " to " + j
        rules[key][0] = row[2]
        rules[key][1] = row[3]
        rules[key][2] = row[4]
        rules[key][3] = row[5]
# Initialise a set of variables to track the metrics calculated.
run_count = 0
# The metrics are stored in a dictionary to account for variable number
# of runs.
metrics = {}
# Specify the working directory, the folder containing the Metronamica
# project file.
working_directory_2000 = ("C:\\Users\\charl\\OneDrive\\Documents\\Geonamica\\"
                          "Metronamica\\" + case_study + "\\")
working_directory_2006 = ("C:\\Users\\charl\\OneDrive\\Documents\\Geonamica\\"
                          "Metronamica\\" + case_study + "_2006\\")
# Specify the project file names.
project_file_2000 = working_directory_2000 + case_study + ".geoproj"
project_file_2006 = working_directory_2006 + case_study + "_2006.geoproj"
# Specify the log files.
log_file_2000 = base_path + "LogSettings.xml"
log_file_2006 = base_path + "LogSettings2006.xml"
# Specify the paths to the simulated output maps.
smap_path_2000 = ("C:\\Users\\charl\\OneDrive\\Documents\\Geonamica\\"
                  "Metronamica\\" + case_study + "\\Log\\Land_use\\"
                  "Land use map_2000-Jan-01 00_00_00.rst")
smap_path_2006 = ("C:\\Users\\charl\\OneDrive\\Documents\\Geonamica\\"
                  "Metronamica\\" + case_study + "_2006\\Log\\Land_use\\"
                  "Land use map_2006-Jan-01 00_00_00.rst")
# Input the rules to the specified project files.
for i in range(0, luc):
    for j in range(0, act):
        key = "from "+luc_names[i]+" to "+luc_names[j+pas]
        fu_elem = j
        lu_elem = i
        y0 = rules[key][0]
        y1 = rules[key][1]
        y2 = rules[key][2]
        xe = rules[key][3]
        set_lp_rule(project_file_2000, fu_elem, lu_elem, y0, y1, y2, xe)
        set_lp_rule(project_file_2006, fu_elem, lu_elem, y0, y1, y2, xe)
# Now generate the simulated output and track the results.
while run_count < max_runs:
    # Perform analysis for the calibration period (1990-2000).
    seed = base_seed + run_count
    set_rand(project_file_2000, seed)
    run_metro(project_file_2000, log_file_2000, working_directory_2000,
              geo_cmd)
    smap = read_map(smap_path_2000)
    it_kappa = kappa(map2, smap, mask)
    it_ksim = ksim(map1, map2, smap, mask)
    it_clu = clumpiness_index(smap, mask, luc)
    key1 = "kappa-2000-it-" + str(seed)
    metrics[key1] = it_kappa
    key2 = "ksim-2000-it-" + str(seed)
    metrics[key2] = it_ksim
    key3 = "clu-2000-it-" + str(seed)
    metrics[key3] = it_clu
    # Perform analysis for validation period (2000-2006).
    set_rand(project_file_2006, seed)
    run_metro(project_file_2006, log_file_2006, working_directory_2006,
              geo_cmd)
    smap = read_map(smap_path_2006)
    it_kappa = kappa(map3, smap, mask)
    it_ksim = ksim(map2, map3, smap, mask)
    it_clu = clumpiness_index(smap, mask, luc)
    key1 = "kappa-2006-it-" + str(seed)
    metrics[key1] = it_kappa
    key2 = "ksim-2006-it-" + str(seed)
    metrics[key2] = it_ksim
    key3 = "clu-2006-it-" + str(seed)
    metrics[key3] = it_clu
    # Add 1 to iterator
    run_count = run_count + 1

# Write the (metrics) to a .csv file

metrics_output_file = output_path + case_study + "\\final_rule_metrics.csv"
store = [0]*(3 + luc)
with open (metrics_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Key", "Kappa", "KSIM"] + luc_names
    writer.writerow(values)
    for x in range(0, max_runs):
        seed = base_seed + x
        key1 = "kappa-2000-it-" + str(seed)
        key2 = "ksim-2000-it-" + str(seed)
        key3 = "clu-2000-it-" + str(seed)
        store[0] = "2000-" + str(x)
        store[1] = metrics[key1]
        store[2] = metrics[key2]
        for b in range(0, luc):
            store[3 + b] = metrics[key3][b]
        writer.writerow(store)
    for x in range(0, max_runs):
        seed = base_seed + x
        key1 = "kappa-2006-it-" + str(seed)
        key2 = "ksim-2006-it-" + str(seed)
        key3 = "clu-2006-it-" + str(seed)
        store[0] = "2006-" + str(x)
        store[1] = metrics[key1]
        store[2] = metrics[key2]
        for b in range(0, luc):
            store[3 + b] = metrics[key3][b]
        writer.writerow(store)
