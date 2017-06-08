"""
This component of the Empirical Neighbourhood Calibration method evaluates the
inertia values and fine-tunes these values to minimise the error between the
simulated and observed enrichment factor values
"""

import numpy as np
import csv
from read_map import read_map
from multi_sims_ef import multi_sims_ef
from considered_distances import considered_distances
from set_NR import set_lp_rule
from log_scale_ef import log_scale_ef

# Specify the base path to the directory containing the empirical
# neighbourhood calibration tool-pack.
base_path = "C:\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Budapest"
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

# Specify the fine-tuning calibration method parameters.
max_runs = 10
base_seed = 1000
# Specify the minimum influence value for inertia points, and the adjustment value.
min_inertia = 200
adj_value = 50


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

# Determine the distances that will be analysed,
# use module: considered_distances.

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
# to store the rules, and input from the initial rule files generated.
rules = {}
att_rules_file = output_path + case_study + "\\Rules\\att_rules.txt"
initial_rules_file = output_path + case_study + "\\Rules\\Initial_rules.csv"
# Initialise an array to track the inertia values
inertia_values = [0]*act
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
        if i == j + pas:
            inertia_values[j] = y0

# Amend the inertia values to generate the conversions present in the data
i_c = np.zeros(shape=(luc, act))
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        if i == j + pas:
            pass
        elif rules[key][0] == str(0):
            pass
        else:
            i_c[i, j] = 1
# Initialise an array, types of data conversions (toc), to count the types of
# meaningful transitions from data analysis.
todc = [0]*act
for i in range(0, act):
    todc[i] = sum(i_c[i + pas, :])

# Begin iterative testing of the Inertia values to simulate the types of
# conversions present in the data.
# Start by calculating the enrichment factor for multiple simulations.
sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                        log_file, working_directory, geo_cmd, luc,
                        max_distance, cdl, cd, N, omap, mask, rows, cols)
# Log scale the simulated enrichment factor values
log_sims_ef = log_scale_ef(sims_ef, 10, luc, act, pas, max_distance)
# Extract the logged enrichment factor values for conversions.
log_sims_ef_0 = log_sims_ef[0, :, :]
# Initialise an array, types of simulated conversions (tosc), to count the types
# of transitions from simulated analysis.
tosc = [0]*act
for i in range(0, act):
    for j in range(0, act):
        if log_sims_ef_0[i + pas, j] != -9999:
            tosc[i] = tosc[i] + 1
# Evaluate the types of conversions that are not occuring.
non_conversions = [0]*act
for i in range(0, act):
    non_conversions[i] = tosc[i] - todc[i]

# Iteratively adjust inertia values to generate types of conversions.
# A dummy variable is initialised to evaluate convergence of the simulation
change = 1
while change > 0:
    # Reset the convervence variable
    change = 0
    # Evaluate the simulated output, and adjust the inertia value for classes
    # where an insufficient number of the types of transitions are occuring.
    for i in range(0, act):
        if non_conversions[i] < 0:
            if float(inertia_values[i]) > float(min_inertia):
                coi = i
                key = "from " + luc_names[i + pas] + " to " + luc_names[i + pas]
                new_y0 = float(rules[key][0]) - adj_value
                inertia_values[i] = new_y0
                rules[key][0] = new_y0
                # Set the new rule value
                fu_elem = coi
                lu_elem = coi + pas
                y0 = rules[key][0]
                y1 = rules[key][1]
                y2 = rules[key][2]
                xe = rules[key][3]
                set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)
                change = 1
    # Re-evaluate the simulated output if an inertia value was adjusted.
    if change == 1:
        sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                                log_file, working_directory, geo_cmd, luc,
                                max_distance, cdl, cd, N, omap, mask, rows, cols)
        log_sims_ef = log_scale_ef(sims_ef, 10, luc, act, pas, max_distance)
        log_sims_ef_0 = log_sims_ef[0, :, :]
        tosc = [0] * act
        for i in range(0, act):
            for j in range(0, act):
                if log_sims_ef_0[i + pas, j] <> -9999:
                    tosc[i] = tosc[i] + 1
        non_conversions = [0] * act
        for i in range(0, act):
            non_conversions[i] = tosc[i] - todc[i]
# Save the generate rules
inertia_adjusted_rules_file = (output_path + case_study +
                               "\\Rules\\inertia_adjusted_rules.csv")
store = [0]*6
with open(inertia_adjusted_rules_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    # Save and write a header line to the file
    values = ["from", "to", "I@x=0", "I@x=1", "I@x=2", "x@I=0"]
    writer.writerow(values)
    # Now write the neighbourhood rules in the form from ... to ...
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            store[0] = luc_names[i]
            store[1] = luc_names[j + pas]
            store[2] = rules[key][0]
            store[3] = rules[key][1]
            store[4] = rules[key][2]
            store[5] = rules[key][3]
            writer.writerow(store)
