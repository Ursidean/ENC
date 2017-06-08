"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine tuning of the neighbourhood rules to improve the agreement with the
observed enrichment factor values.
"""

import numpy as np
from read_map import read_map
from considered_distances import considered_distances
import csv
from set_NR import set_lp_rule
from enrichment_factor import ef
from log_scale_ef import log_scale_ef
from multi_sims_ef import multi_sims_ef
import math
from numpy import unravel_index
import time


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Berlin"
# Set the paths to the relevant directories
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
working_directory = ("C:\\Users\\a1210607\\Geonamica\\Metronamica\\" +
                     case_study + "\\")
# Set the paths to the relevant files & executables
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
smap_path = working_directory + ("Log\\Land_use\\"
                                 "Land use map_2000-Jan-01 00_00_00.rst")
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
project_file = working_directory + case_study + ".geoproj"
log_file = base_path + "LogSettings.xml"

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
max_runs = 3
base_seed = 1000
gss_tol = 1.0
# Set the base that enrichment factor values are log-scaled to. Standard is 10.
log_base = 10

# Read in the map for the data at time slice 1.
omap = read_map(map1_path)
# Read in the map for the data at time slice 2.
amap = read_map(map2_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
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

# Initialise an array to store the requisite attraction rules and a dictionary
# to store the rules, and input from the inertia adjusted rule file generated.
rules = {}
att_rules_file = output_path + case_study + "\\Rules\\att_rules.txt"
initial_rules_file = (output_path + case_study +
                      "\\Rules\\inertia_adjusted_rules.csv")
# Input the attraction rules array
att_rules = np.loadtxt(att_rules_file)
# Input the initial rules dictionary
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        rules[key] = [0, 0, 0, max_distance]
with open(initial_rules_file, 'rb') as f:
    # Skip the header line in the csv file
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
# Set the neighbourhood rules in the Metronamica geoproject file.
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        fu_elem = j
        lu_elem = i
        y0 = rules[key][0]
        y1 = rules[key][1]
        y2 = rules[key][2]
        xe = rules[key][3]
        set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)

# Golden Section Search (GSS) variable initialisation.
# Determine the enrichment of the data, and log scale the values.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
log_data_ef = log_scale_ef(data_ef, log_base, luc, act, pas, max_distance)
# Run an initial simulation to generate the initial simulated enrichment
# factor values, and log-scale the values
sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                        log_file, working_directory, geo_cmd, luc,
                        max_distance, cdl, cd, N, omap, mask, rows, cols)
log_sims_ef = log_scale_ef(sims_ef, log_base, luc, act, pas, max_distance)
# Initialise a set of empty lists:
# 1. ssr_sum_log is used to track the error convergance of the GSS method.
# 2. rule_tracker is used to track the rule being adjusted.
ssr_sum_log = []
rule_tracker = []
# Initialise a variable to track the best rules obtained.
best_rules = rules
# Generate the initial Sum of Square Residuals (SSR) values

