"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine tuning of the neighbourhood rules to improve the agreement with the
observed enrichment factor values.
"""

import math
import numpy as np
from read_map import read_map
from considered_distances import considered_distances
import csv
from set_NR import set_lp_rule
from enrichment_factor import ef
from log_scale_ef import log_scale_ef
from multi_sims_ef import multi_sims_ef
import time
from numpy import unravel_index
from ef_ssr_calc import ef_ssr
from kappa import ksim


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
# maximum number of simulation runs, the golden section search tolerance, and
# an error value for when the sum of square residuals can't be evaluated
# (default value of 5.0).
max_runs = 3
base_seed = 1000
gss_tol = 2.5
no_ssr_eval_error = 5.0

# Specify the golden ratio (approx.).
gr = (math.sqrt(5) + 1)/2
# Specify the maximum and minimum influence values. Defaults are 0 and 100
# for both conversion points and influence at a distance 1.
max_conversion = 100
min_y1_influence = 0
max_y1_influence = 100
# Set the base that enrichment factor values are log-scaled to. Default is 10.
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
# factor values, and log-scale the values.
sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                        log_file, working_directory, geo_cmd, luc,
                        max_distance, cdl, cd, N, omap, mask, rows, cols)
log_sims_ef = log_scale_ef(sims_ef, log_base, luc, act, pas, max_distance)
# Initialise a set of empty lists:
# 1. ssr_sum_log is used to track the error convergence of the GSS method.
# 2. rule_tracker is used to track the rule being adjusted.
# 3. ksim_log tracks the kappa simulation values for each iteration.
ssr_sum_log = []
rule_tracker = []
ksim_log = []
# Initialise a variable to track the best rules obtained.
best_rules = rules
# Generate the initial Sum of Square Residuals (SSR) values

dummy = ef_ssr(log_data_ef, log_sims_ef, luc, act, no_ssr_eval_error)
ssr_matrix = dummy[0]
ssr_sum = dummy[1]
ssr_sum_log.append(ssr_sum)
smap = read_map(smap_path)
initial_ksim = ksim(omap, amap, smap, mask)
ksim_log.append(initial_ksim)

# Set the previous iterations SSR value to initialise the calibration.
ssr_sum_old = ssr_sum + 1
# Initialise a set of variables to prevent successive rule tuning
lu_elem_old = luc - 1
fu_elem_old = act - 1
# Begin a timer to track the duration of the calibration method
start_ft = time.time()

# Golden Section Search iterative adjustment
while ssr_sum < ssr_sum_old:
    # Set the old SSR to the previous iteration value
    ssr_sum_old = ssr_sum
    # Set the previous rule SSR to 0 to prevent consecutive adjustment.
    ssr_matrix[lu_elem_old, fu_elem_old] = 0
    # Determine which rule to adjust. Start by setting non-evaluated
    # values to error code so these are not adjusted.
    for i in range(0, luc):
        for j in range(0, act):
            if att_rules[i, j] != 1:
                ssr_matrix[i, j] = -9999
            elif ssr_matrix[i, j] == no_ssr_eval_error:
                ssr_matrix[i, j] = -9999
    # Now determine the index of the rule to adjust.
    rule_index = unravel_index(ssr_matrix.argmax(), ssr_matrix.shape)
    lu_elem = rule_index[0]
    fu_elem = rule_index[1]
    rule_name = ("from " + luc_names[lu_elem] +
                 " to " + luc_names[fu_elem + pas])
    print "Calibrating rule: " + rule_name
    # Find the influence values for the specified rule
    adjust_rule_values = rules[rule_name]
    # Fix the ratio between the influence values at distance 1 & 2.
    ratio_21 = float(adjust_rule_values[2]) / float(adjust_rule_values[1])
    # Stage 1. Set to extreme values and initialise the GSS.
    # Set fixed elements for the neighbourhood rule adjustment.
    y0 = adjust_rule_values[0]
    xe = adjust_rule_values[3]
    # Specify the low influence values for the required distances.
    y1_low = min_y1_influence
    y2_low = min_y1_influence*ratio_21
    # Specify the high influence values for the required distances.
    y1_high = max_y1_influence
    y2_high = max_y1_influence * ratio_21
    # Initialise the bracketing influence values, low (a) and high (b).
    y1_a = y1_low
    y2_a = y2_low
    y1_b = y1_high
    y2_b = y2_high
    # Stage 2. Apply the Golden Section Search to converge on the
    # best influence values
    while (abs(y1_a - y1_b) > gss_tol):
        y1_c = y1_b - (y1_b - y1_a)/gr
        y1_d = y1_a + (y1_b - y1_a)/gr
        y2_c = y1_c*ratio_21
        y2_d = y1_d*ratio_21
        # Set the neighbourhood rule to the lower interval values.
        set_lp_rule(project_file, fu_elem, lu_elem, y0, y1_c, y2_c, xe)
        # Calculate the simulated enrichment factors, log scale.
        c_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path,
                                  project_file, log_file,
                                  working_directory, geo_cmd, luc,
                                  max_distance, cdl, cd, N, omap,
                                  mask, rows, cols)
        log_c_sims_ef = log_scale_ef(c_sims_ef, log_base, luc, act, pas,
                                     max_distance)
        # Calculate the sum of square residuals for the lower interval.
        ssr_c = ef_ssr(log_data_ef, log_c_sims_ef, luc, act,
                       no_ssr_eval_error)
        f_c = ssr_c[1]
        # Set the neighbourhood rule to the higher interval values.
        set_lp_rule(project_file, fu_elem, lu_elem, y0, y1_d, y2_d, xe)
        # Calculate the simulated enrichment factors, log scale.
        d_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path,
                                  project_file, log_file,
                                  working_directory, geo_cmd, luc,
                                  max_distance, cdl, cd, N, omap,
                                  mask, rows, cols)
        log_d_sims_ef = log_scale_ef(d_sims_ef, log_base, luc, act, pas,
                                     max_distance)
        # Calculate the sum of square residuals for the higher interval.
        ssr_d = ef_ssr(log_data_ef, log_d_sims_ef, luc, act,
                       no_ssr_eval_error)
        f_d = ssr_d[1]
        # Print output
        print ("Current influence bracket: " + str(y1_c) + " | " +
               str(y1_d))
        print ("SSR influence bracket: " + str(f_c) + " | " + str(f_d))
        # Evaluate the result
        if f_c < f_d:
            y1_b = y1_d
            y2_b = y2_d
        else:
            y1_a = y1_c
            y2_a = y2_c
    # Process the final evlauation, and set to final rule.
    if f_c <= f_d:
        rules[rule_name] = [y0, y1_c, y2_c, xe]
    else:
        rules[rule_name] = [y0, y1_d, y2_d, xe]
    # Set the final values.
    set_lp_rule(project_file, fu_elem, lu_elem, rules[rule_name][0],
                rules[rule_name][1], rules[rule_name][2],
                rules[rule_name][3])
    # Re-calculate the enrichment factors, log scale.
    sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                            log_file, working_directory, geo_cmd, luc,
                            max_distance, cdl, cd, N, omap, mask, rows, cols)
    log_sims_ef = log_scale_ef(sims_ef, log_base, luc, act, pas, max_distance)
    ssr_output = ef_ssr(log_data_ef, log_sims_ef, luc, act,
                        no_ssr_eval_error)
    # Set the new SSR value, log the values.
    ssr_matrix = ssr_output[0]
    ssr_sum = ssr_output[1]
    ssr_sum_log.append(ssr_sum)
    # Evaluate the iteration Kappa Simulation
    #smap = read_map(smap_path)
    #it_ksim = ksim(omap, amap, smap, mask)

    # Evaluate and track the rule adjusted.
    rule_tracker.append(rule_name)
    rule_tracker.append(rules[rule_name])
    # Track the function and land-use element to prevent consecutive
    # adjustment.
    lu_elem_old = lu_elem
    fu_elem_old = fu_elem
    # Track the best rules
    if ssr_sum < ssr_sum_old:
        best_rules = rules

end_ft = time.time()

# Store the output generated from the fine-tuning procedure.
# Print the time taken.
print "Time taken:" + str(end_ft - start_ft)
# Log the final rules set.
final_rules_file_path = output_path + case_study + "\\Rules\\final_rules.csv"
store = [0]*6
with open(final_rules_file_path, "wb") as csv_file:
    writer = csv.writer(csv_file)
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
# Log the sum of square residuals
ssr_log_file = output_path + case_study + "\\ssr_log.csv"
store = [0]*2
with open(ssr_log_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Iteration", "SSR sum", "KSIM"]
    writer.writerow(values)
    for i in range(0, len(ssr_sum_log)):
        store[0] = i
        store[1] = ssr_sum_log[i]
        #store[2] = ksim_log[i]
        writer.writerow(store)
# Log the specific adjustments.
rule_adjustment_file = (output_path + case_study +
                        "\\Rules\\ft_rule_adjusted.txt")
thefile = open(rule_adjustment_file, "w")
for item in rule_tracker:
    thefile.write("%s\n" % item)
thefile.close()
# Fine tuning completed!
