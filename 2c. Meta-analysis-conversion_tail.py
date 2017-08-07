"""
This component of the Empirical Neighbourhood Calibration method conducts the
enumeration analysis of the meta-parameters used to initialise the
calibration method.
"""

import numpy as np
from considered_distances import considered_distances
from read_map import read_map
from enrichment_factor import ef
from log_scale_ef import log_scale_ef
from contingency_table import contingency_table
from set_NR import set_lp_rule
from set_rand import set_rand
from run_metro import run_metro
from kappa import ksim
from area_weighted_clu import area_weighted_clu_error
import csv

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Specify the case study.
case_study = "Budapest"
# Set the paths to the directories and relevant data
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the working directory (location of project file) and project file.
working_directory = ("C:\\Geonamica\\Metronamica\\"
                              + case_study + "\\")
project_file = working_directory + case_study + ".geoproj"
# Specify the paths to the simulated output maps.
smap_path = ("C:\\Geonamica\\Metronamica\\"
                      + case_study + "\\Log\\Land_use\\"
                                     "Land use map_2000-Jan-01 00_00_00.rst")
# Read in the data maps and mask.
omap = read_map(map1_path)
amap = read_map(map2_path)
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]
# Specify the command line version of Geonamica.
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Specify the log file path.
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

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1
# Specify the maximum neighbourhood size distance considered
max_distance = 5
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

# Specify the meta-analysis calibration method parameters: the base random
# seed, and the maximum number of simulation runs.
max_runs = 10
base_seed = 1000

# Specify the bands levels.
# The inertia band levels, for setting inertia values.
high_inertia_band = 0.95
mid_inertia_band = 0.90
# The conversion band levels, for setting conversion values.
high_conversion_band = 0.5
mid_conversion_band = 0.1
low_conversion_band = 0.025
# The enrichment factor band levels, for setting the tail values.
high_ef = 1.0
mid_ef = 0.5

# Specify high, medium and low inertia values.
# The default settings are a high inertia of 1000, med 500, low 250.
high_inertia = 1000.0
mid_inertia = 500.0
low_inertia = 250.0

"""***Change the values here ***"""
# Set default meta-parameters.
# Self-influence tail, specified at a distance of 1
theta_si1 = 0.025
theta_si2 = theta_si1*0.1
# Conversion point.
theta_cp = 0.045
# Interaction tail, specified at a distance of 1
theta_ct1 = 0.005
theta_ct2 = theta_ct1*0.1

# Calculate all relevant neighbourhood rule parameter values.
# Set the conversion parameter values.
high_conversion = high_inertia * theta_cp
mid_conversion = mid_inertia * theta_cp
low_conversion = low_inertia * theta_cp
# Set the distance 1 tail values for self-influence rules.
d1_high_si_value = high_inertia * theta_si1
d1_mid_si_value = mid_inertia * theta_si1
d1_low_si_value = low_inertia * theta_si1
# Set the distance 2 tail values for self-influence rules.
d2_high_si_value = high_inertia * theta_si2
d2_mid_si_value = mid_inertia * theta_si2
d2_low_si_value = low_inertia * theta_si2
# Set the distance 1 tail values for interaction rules.
d1_high_co_value = high_inertia * theta_ct1
d1_mid_co_value = mid_inertia * theta_ct1
d1_low_co_value = low_inertia * theta_ct1
# Set the distance 2 tail value for interaction rules.
d2_high_co_value = high_inertia * theta_ct2
d2_mid_co_value = mid_inertia * theta_ct2
d2_low_co_value = low_inertia * theta_ct2

# Generate the Enrichment Factor and contingency table.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
# Log scale the enrichment factor values.
log_data_ef = log_scale_ef(data_ef, 10, luc, act, pas, max_distance)
# Generate the contingency table using the module, 'contingency_table'
cont_table = contingency_table(omap, amap, mask, luc, rows, cols)
# Evaluate the rates of inertia and conversion
ic_rates = np.zeros(shape=(luc, luc))
for i in range(0, luc):
    for j in range(0, luc):
        if i == j:
            if cont_table[i, luc] > 0:
                ic_rates[i, j] = cont_table[i, j] / cont_table[i, luc]
        else:
            conversions = abs(float(cont_table[j, j]) - float(cont_table[luc, j]))
            if conversions > 0:
                ic_rates[i, j] = float(cont_table[i, j]) / float(conversions)
# Load the attraction rules file
att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
att_rules = np.loadtxt(att_rule_file)
# Input the rules to be analysed. Start by initialising a dictionary for
# storage.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, 5]
# Set the initial neighbourhood rule values for inertia and conversion.
for i in range(0, act):
    for j in range(0, luc):
        # Specify the neighbourhood rule key.
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        # If a self-influence rule, set the inertia value.
        if i + pas == j:
            if cont_table[i + pas, luc] > cont_table[luc, i + pas]:
                rules[key][0] = low_inertia
            else:
                inertia_rate = ic_rates[j, i + pas]
                if inertia_rate > high_inertia_band:
                    rules[key][0] = high_inertia
                elif inertia_rate > mid_inertia_band:
                    rules[key][0] = mid_inertia
                else:
                    rules[key][0] = low_inertia
        # If an interactive rule, set the conversion rule.
        else:
            conversion_rate = ic_rates[j, i + pas]
            if conversion_rate > high_conversion_band:
                rules[key][0] = high_conversion
            elif conversion_rate > mid_conversion_band:
                rules[key][0] = mid_conversion
            elif conversion_rate > low_conversion_band:
                rules[key][0] = low_conversion
# Set the initial neighbourhood rule values for attraction.
for i in range(0, act):
    for j in range(0, luc):
        # Specify the neighbourhood rule key.
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        # If a self-influence rule, set the self-influence attraction values.
        if i + pas == j:
            for c in range(1, 3):
                if c == 1:
                    if att_rules[j, i] == 1:
                        if log_data_ef[c, j, i] > high_ef:
                            rules[key][c] = d1_high_si_value
                        elif log_data_ef[c, j, i] > mid_ef:
                            rules[key][c] = d1_mid_si_value
                        else:
                            rules[key][c] = d1_low_si_value
                elif c == 2:
                    if att_rules[j, i] == 1:
                        if log_data_ef[c, j, i] > high_ef:
                            rules[key][c] = d2_high_si_value
                        elif log_data_ef[c, j, i] > mid_ef:
                            rules[key][c] = d2_mid_si_value
                        else:
                            rules[key][c] = d2_low_si_value
        # If a conversion rule, set the interactive attraction values.
        else:
            if (
                att_rules[j, i] == 1 and log_data_ef[1, j, i] > 0 and
                log_data_ef[2, j, i] > 0
            ):
                for c in range(1, 3):
                    if c == 1:
                        if log_data_ef[c, j, i] > high_ef:
                            rules[key][c] = d1_high_co_value
                        elif log_data_ef[c, j, i] > mid_ef:
                            rules[key][c] = d1_mid_co_value
                        elif log_data_ef[c, j, i] > 0:
                            rules[key][c] = d1_low_co_value
                    elif c == 2:
                        if log_data_ef[c, j, i] > high_ef:
                            rules[key][c] = d2_high_co_value
                        elif log_data_ef[c, j, i] > mid_ef:
                            rules[key][c] = d2_mid_co_value
                        elif log_data_ef[c, j, i] > 0:
                            rules[key][c] = d2_low_co_value
# Set the end-points of each attraction rule
for i in range(0, act):
    for j in range(0, luc):
        if att_rules[j, i] == 0:
            pass
        else:
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # Iterate through to find end point
            for c in range(2, 5):
                if att_rules[j, i] == 1 and log_data_ef[c, j, i] > 0:
                    rules[key][3] = c + 1

# Input the rules into the model.
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

# Enumerate possible ratio values.
# First, initialise a dictionary to store metrics, this is done to account
# for variable size.
metrics = {}
"""***Change the bounds here ***"""
# Start with the self-influence tail values. First, specify the range that will be tested.
varied_parameter = "theta_ct"
low_bound = 0
high_bound = 1
theta_ct_step = 0.001
for a in range(low_bound, high_bound, 1):
    theta_ct1 = float(a)*theta_ct_step
    theta_ct2 = theta_ct1*0.1
    # Set the distance 1 tail values for interaction rules.
    d1_high_co_value = high_inertia * theta_ct1
    d1_mid_co_value = mid_inertia * theta_ct1
    d1_low_co_value = low_inertia * theta_ct1
    # Set the distance 2 tail value for interaction rules.
    d2_high_co_value = high_inertia * theta_ct2
    d2_mid_co_value = mid_inertia * theta_ct2
    d2_low_co_value = low_inertia * theta_ct2
    # Input the set values into the model.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, skip.
            if i + pas == j:
                pass
            # If an interactive rule, set the conversion point value.
            else:
                if (
                    att_rules[j, i] == 1 and log_data_ef[1, j, i] > 0 and
                    log_data_ef[2, j, i] > 0
                ):
                    for c in range(1, 3):
                        if c == 1:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d1_high_co_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d1_mid_co_value
                            elif log_data_ef[c, j, i] > 0:
                                rules[key][c] = d1_low_co_value
                        elif c == 2:
                            if log_data_ef[c, j, i] > high_ef:
                                rules[key][c] = d2_high_co_value
                            elif log_data_ef[c, j, i] > mid_ef:
                                rules[key][c] = d2_mid_co_value
                            elif log_data_ef[c, j, i] > 0:
                                rules[key][c] = d2_low_co_value
    # Input the rules to the specified project files.
    # This is only done for the conversion rules.
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            fu_elem = j
            lu_elem = i
            if i != j + pas:
                y0 = rules[key][0]
                y1 = rules[key][1]
                y2 = rules[key][2]
                xe = rules[key][3]
                set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)
    # Generate the simulated output and track the results.
    run_count = 0
    ksim_log = [0]*max_runs
    clu_log = [0]*max_runs
    while run_count < max_runs:
        # Generate seed, run model to generate output.
        rseed = base_seed + run_count
        set_rand(project_file, rseed)
        run_metro(project_file, log_file, working_directory,
                  geo_cmd)
        # Read in the map.
        smap = read_map(smap_path)
        # Calculate the corresponding metrics
        ksim_log[run_count] = ksim(omap, amap, smap, mask)
        clu_log[run_count] = area_weighted_clu_error(amap, smap, mask, luc, pas, act, luc_count)
        # Add one to iterator to avoid infinite loop!
        run_count = run_count + 1
    # Log metric properties corresponding to the output
    # First for Kappa Simulation
    key = "ksim-" + str(a)
    metrics[key] = [0]*3
    metrics[key][0] = sum(ksim_log)/len(ksim_log)
    metrics[key][1] = (max(ksim_log) - min(ksim_log))/2
    # Next for Absolute average clumpiness error.
    key = "clu_error-" + str(a)
    metrics[key] = [0] * 3
    metrics[key][0] = sum(clu_log) / len(clu_log)
    metrics[key][1] = (max(clu_log) - min(clu_log)) / 2
# Organise and evaluate output.

# Write the output to a .csv file.
metrics_output_file = output_path + case_study + "\\theta_ct_meta_analysis_output.csv"
store = [0]*5
with open (metrics_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["theta_ct value", "KSIM avg.", "KSIM_err", "CLU avg", "CLU_err"]
    writer.writerow(values)
    for a in range(low_bound, high_bound, 1):
        store[0] = float(a)*theta_ct_step
        store[1] = metrics[("ksim-" + str(a))][0]
        store[2] = metrics[("ksim-" + str(a))][1]
        store[3] = metrics[("clu_error-" + str(a))][0]
        store[4] = metrics[("clu_error-" + str(a))][1]
        writer.writerow(store)

# Finished!