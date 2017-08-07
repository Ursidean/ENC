"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine tuning of the neighbourhood rules to improve the agreement with the
observed enrichment factor values.
"""

import math
import numpy as np
from read_map import read_map
from considered_distances import considered_distances
from contingency_table import contingency_table
from enrichment_factor import ef
from log_scale_ef import log_scale_ef
from set_NR import set_lp_rule
import operator
from set_rand import set_rand
from run_metro import run_metro
from kappa import ksim
from area_weighted_clu import area_weighted_clu_error
from multi_sims_ef import multi_sims_ef
import time
from numpy import unravel_index
from ef_ssr_calc import ef_ssr
import csv

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Budapest"
# Set the paths to the relevant directories
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
working_directory = ("C:\\Geonamica\\Metronamica\\" +
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

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1

# Part 1. Initial rule generation.
# Generate an initial rule set to start the calibration.
# Specify the band settings for determining the initial neighbourhood rules.
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
# Specify the high, medium and low inertia point values.
high_inertia = 1000.0
mid_inertia = 500.0
low_inertia = 250.0
# Specify the meta-parameter values.
theta_i1 = 0.025
theta_i2 = theta_i1*0.1
theta_cp = 0.045
theta_c1 = 0.003
theta_c2 = theta_c1*0.1
# Calculate the relevant neighbourhood rule parameter values.
# Set the distance 1 tail values for self-influence rules.
d1_high_si_value = high_inertia * theta_i1
d1_mid_si_value = mid_inertia * theta_i1
d1_low_si_value = low_inertia * theta_i1
# Set the distance 2 tail values for self-influence rules.
d2_high_si_value = high_inertia * theta_i2
d2_mid_si_value = mid_inertia * theta_i2
d2_low_si_value = low_inertia * theta_i2
# Set the conversion parameter values.
high_conversion = high_inertia * theta_cp
mid_conversion = mid_inertia * theta_cp
low_conversion = low_inertia * theta_cp
# Set the distance 1 tail values for interaction rules.
d1_high_co_value = high_inertia * theta_c1
d1_mid_co_value = mid_inertia * theta_c1
d1_low_co_value = low_inertia * theta_c1
# Set the distance 2 tail value for interaction rules.
d2_high_co_value = high_inertia * theta_c2
d2_mid_co_value = mid_inertia * theta_c2
d2_low_co_value = low_inertia * theta_c2
# Input the initial rules into the model.
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
            if (att_rules[j, i] == 1 and log_data_ef[1, j, i] > 0 and
                        log_data_ef[2, j, i] > 0):
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

# Part 2. Fine tuning of parameters.
# Successively fine tune the individual parameters of the model.
# Specify the golden ratio (approx).
gr = (math.sqrt(5) + 1)/2
# Specify the parameters of the golden section search.
max_runs = 10
base_seed = 1000
# Specify the maximum and minimum bounding values for different neighbourhood
# components, and the golden section search tolerance:
# Inertia point
max_ip = 1000
min_ip = 250
gss_ip_tol = (max_ip - min_ip)/100
# The self-influence tail influence value at distance 1.
max_si = 100
min_si = 0
gss_si_tol = (max_si - min_si)/100
# The conversion point
max_cp = 100
min_cp = 0
gss_cp_tol = (max_cp - min_cp)/100
# The interaction tail at distance 1.
max_ct = 100
min_ct = 0
gss_ct_tol = (max_ct - min_ct)/100
# Initialise a variable to track the best rules obtained.
best_rules = rules
# Initialise two lists to track the rule adjusted by name.
rule_from_tracker = []
rule_to_tracker = []
# Initialise a list to track the point that is fine tuned,
# and the final value taken.
pt_tracker = []
pt_value = []
# Initialise an array to track the values of inertia points generated.
ip_value_tracker = [0]*act

# Set the order to calibrate the variables. Start by specifying the
# socio-economic group of each land-use class. This will be used to
# set the order.
socio_eco_group = ["Other", "Other", "Work", "Work", "Work", "Residential",
                   "Work", "Recreation", "Other", "Other", "Other", "Work",
                   "Work", "Other", "Other", "Other", ]
# Initialise a pair of lists to track the area-weighted absoluate average
# clumpiness error and the Kappa Simulation value.
clu_log = []
ksim_log= []
# Determine the present Kappa Simulation.
current_ksim = [0]*max_runs
current_ce = [0]*max_runs
run_count = 0
while run_count < max_runs:
    # Set the random seed
    rseed = base_seed + run_count
    set_rand(project_file, rseed)
    # Run Metronamica to generate output
    run_metro(project_file, log_file, working_directory,
              geo_cmd)
    # Read in output
    smap = read_map(smap_path)
    # Determine KSIM
    current_ksim[run_count] = ksim(omap, amap, smap, mask)
    current_ce[run_count] = area_weighted_clu_error(amap,smap, mask, luc, pas,
                                                    act, luc_count)
    # Add one to iterator to prevent infinite loop!
    run_count = run_count + 1
ksim_avg = sum(current_ksim)/len(current_ksim)
clu_avg = sum(current_ce)/len(current_ce)
# Append the values to the list.
#ksim_log.append(ksim_avg)
#clu_log.append(clu_avg)




# 2a. Fine tuning of inertia points. ---
print "Fine tuning the inertia points"
# Initialise a list for self-influence rule order tuning.
ip_order = np.array([0]*act)
# Evaluate the socio-economic groups to determine the order for
# inertia points and self-influence rules.
# Allocate a numerical value.
for i in range(0, act):
    if socio_eco_group[i + pas] == "Residential":
        ip_order[i] = 4
    elif socio_eco_group[i + pas] == "Work":
        ip_order[i] = 3
    elif socio_eco_group[i + pas] == "Recreation":
        ip_order[i] = 2
    else:
        ip_order[i] = 1
# Rank and determine order.
temp = ip_order.argsort()
ranks = np.empty(len(ip_order), int)
ranks[temp] = np.arange(len(ip_order))
ip_order = ranks + 1
# Apply the Golden Section Search to Inertia points.
while sum(ip_order) > 0:
    # Identify max value, corresponding to parameter that will be
    # tested.
    par_index, par_value  = max(enumerate(ip_order), key=operator.itemgetter(1))
    # Set the corresponding value in the list to 0.
    ip_order[par_index] = 0
    # Specify the function and land-use element index values for
    # appending the geoproject file.
    fu_elem = par_index
    lu_elem = par_index + 1
    rule_name = ("from " + luc_names[lu_elem] +
                 " to " + luc_names[fu_elem + pas])
    print "Calibrating rule: " + rule_name
    # Track the rule being adjusted
    rule_from_tracker.append(lu_elem)
    rule_to_tracker.append(fu_elem)
    pt_tracker.append(0)
    # Find the influence values for the specified rule
    adjust_rule_values = rules[rule_name]
    # Fix the requisite parameter values.
    old_y0 = adjust_rule_values[0]
    y1 = adjust_rule_values[1]
    y2 = adjust_rule_values[2]
    xe = adjust_rule_values[3]
    # Apply the Golden Section Search to optimise the output.
    a = min_ip
    b = max_ip
    while (abs(a - b) > gss_ip_tol):
        # Calculate the interval values.
        c = b - (b - a)/gr
        d = a + (b - a)/gr
        # Set the neighbourhood rule to the lower interval values.
        set_lp_rule(project_file, fu_elem, lu_elem, c, y1, y2, xe)
        # Run the model, evaluate performance.
        ksim_c = [0]*max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_c[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_c = sum(ksim_c)/len(ksim_c)
        # Now set the neighbourhood rule to the higher interval value.
        set_lp_rule(project_file, fu_elem, lu_elem, d, y1, y2, xe)
        # Run the model, evaluate performance.
        ksim_d = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_d[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_d = sum(ksim_d)/len(ksim_d)
        # Add step here to break out if not ksim not changing.
        if ksim_avg == ksim_c and ksim_avg == ksim_d:
            break
        else:
            # Evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            if ksim_c > ksim_d:
                b = d
            else:
                a = c
    # Set the final value, input into the model.
    if ksim_c == ksim_d and ksim_avg == ksim_c and ksim_avg == ksim_d:
        rules[rule_name][0] = old_y0
        set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
        ksim_log.append(ksim_avg)
        clu_log.append(clu_avg)
        pt_value.append(old_y0)
    else:
        final_value = (a + b)/2
        set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
        # Evaluate the final value.
        ksim_final = [0]*max_runs
        clu_final = [0]*max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_final[run_count] = ksim(omap, amap, smap, mask)
            clu_final[run_count] = area_weighted_clu_error(amap,smap, mask,
                                                           luc, pas, act,
                                                           luc_count)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine whether to change rule or not based on result.
        ksim_avg_new = sum(ksim_final)/len(ksim_final)
        clu_avg_new = sum(clu_final)/len(clu_final)
        if ksim_avg_new < ksim_avg:
            # Conditional. If worse performance, reset.
            rules[rule_name][0] = old_y0
            set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)
            pt_value.append(old_y0)
        else:
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            ksim_avg = ksim_avg_new
            clu_avg = clu_avg_new
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)





# 2b. Fine tuning of self-influence tails. ---
print "Fine tuning the self-influence tails"
# Initialise a list for self-influence rule order tuning.
si_order = np.array([0]*act)
# Evaluate the socio-economic groups to determine the order for
# inertia points and self-influence rules.
# Allocate a numerical value.
for i in range(0, act):
    if socio_eco_group[i + pas] == "Residential":
        si_order[i] = 4
    elif socio_eco_group[i + pas] == "Work":
        si_order[i] = 3
    elif socio_eco_group[i + pas] == "Recreation":
        si_order[i] = 2
    else:
        si_order[i] = 1
# Rank and determine order.
temp = si_order.argsort()
ranks = np.empty(len(si_order), int)
ranks[temp] = np.arange(len(si_order))
si_order = ranks + 1
# Apply the Golden Section Search to self-influence tails.
while sum(si_order) > 0:
    # Identify max value, corresponding to parameter that will be
    # tested.
    par_index, par_value = max(enumerate(si_order), key=operator.itemgetter(1))
    # Set the corresponding value in the list to 0.
    si_order[par_index] = 0
    # Specify the function and land-use element index values for
    # appending the geoproject file.
    fu_elem = par_index
    lu_elem = par_index + 1
    rule_name = ("from " + luc_names[lu_elem] +
                 " to " + luc_names[fu_elem + pas])
    print "Calibrating rule: " + rule_name
    # Track the rule being adjusted
    rule_from_tracker.append(lu_elem)
    rule_to_tracker.append(fu_elem)
    pt_tracker.append(1)
    # Find the influence values for the specified rule
    adjust_rule_values = rules[rule_name]
    # Fix the requisite parameter values.
    y0 = adjust_rule_values[0]
    old_y1 = adjust_rule_values[1]
    old_y2 = adjust_rule_values[2]
    xe = adjust_rule_values[3]
    # Apply the Golden Section Search to optimise the output.
    a = min_si
    b = max_si
    while (abs(a - b) > gss_si_tol):
        # Calculate the interval values.
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        # Set the neighbourhood rule to the lower interval values.
        c2 = 0.1*c
        set_lp_rule(project_file, fu_elem, lu_elem, y0, c, c2, xe)
        # Run the model, evaluate performance.
        ksim_c = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_c[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_c = sum(ksim_c) / len(ksim_c)
        # Now set the neighbourhood rule to the higher interval value.
        d2 = 0.1 * d
        set_lp_rule(project_file, fu_elem, lu_elem, y0, d, d2, xe)
        # Run the model, evaluate performance.
        ksim_d = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_d[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_d = sum(ksim_d) / len(ksim_d)
        # Add step here to break out if not ksim not changing.
        if ksim_avg == ksim_c and ksim_avg == ksim_d:
            break
        else:
            # Evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            if ksim_c > ksim_d:
                b = d
            else:
                a = c
    # Set the final value, input into the model.
    if ksim_avg == ksim_c and ksim_avg == ksim_d:
        rules[rule_name][1] = old_y1
        rules[rule_name][2] = old_y2
        set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1,
                    old_y2, xe)
        ksim_log.append(ksim_avg)
        clu_log.append(clu_avg)
        pt_value.append(old_y1)
    else:
        final_value = (a + b) / 2
        final_value_2 = final_value * 0.1
        rules[rule_name][1] = final_value
        rules[rule_name][2] = final_value_2
        set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                    final_value_2, xe)
        # Evaluate the final value.
        ksim_final = [0] * max_runs
        clu_final = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM & clumpiness error.
            ksim_final[run_count] = ksim(omap, amap, smap, mask)
            clu_final[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                           luc, pas, act,
                                                           luc_count)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine whether to change rule or not based on result.
        ksim_avg_new = sum(ksim_final)/len(ksim_final)
        clu_avg_new = sum(clu_final)/len(clu_final)
        if ksim_avg_new < ksim_avg:
            # Conditional. If worse performance, reset.
            rules[rule_name][1] = old_y1
            rules[rule_name][2] = old_y2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1,
                        old_y2, xe)
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)
            pt_value.append(old_y1)
        else:
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            ksim_avg = ksim_avg_new
            clu_avg = clu_avg_new
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)
            pt_value.append(final_value)





# 2c. Fine tuning of conversion points. ---
print "Fine tuning the conversion points"
temp_cp_order = np.zeros(shape=(luc,act))
# Evaluate the socio-economic groups to determine the order for
# conversion points. First, allocate a numerical value depending on
# class converted to.
for i in range(0, luc):
    for j in range(0, act):
        if i == j + pas:
            # Skip if an inertia point.
            pass
        elif ic_rates[i, j + pas] < low_conversion_band:
            # Skip if no conversion point in model
            pass
        else:
            # Evaluate the conversions, assign rank based on class that is
            # being converted to.
            if socio_eco_group[j + pas] == "Residential":
                temp_cp_order[i, j] = 4
            elif socio_eco_group[j + pas] == "Work":
                temp_cp_order[i, j] = 3
            elif socio_eco_group[j + pas] == "Recreation":
                temp_cp_order[i, j] = 2
            else:
                temp_cp_order[i, j] = 1
# Rank and determine order. Start by initialising an array to track order.
cp_order = np.zeros(shape=(luc, act))
# Iterate through and assign rank. Start at 1.
ranking = 1
iterator = 1
while iterator < 5:
    # Iterate through the ranks.
    for i in range(0, luc):
        for j in range(0, act):
            if temp_cp_order[i, j] == iterator:
                # If the equivalent value is found, assign current
                # ranking value.
                cp_order[i, j] = ranking
                # Add one to the ranking value.
                ranking = ranking + 1
    # Add one to the iterator.
    iterator = iterator + 1
# Apply the Golden-Section Search to optimise the output.
while np.sum(cp_order > 0):
    # Identify max value, corresponding to parameter that will be
    # tested.
    par_index = unravel_index(cp_order.argmax(), cp_order.shape)
    # Set the corresponding value in the list to 0.
    cp_order[par_index[0], par_index[1]] = 0
    # Specify the function and land-use element index values for
    # appending the geoproject file.
    lu_elem = par_index[0]
    fu_elem = par_index[1]
    rule_name = ("from " + luc_names[lu_elem] +
                 " to " + luc_names[fu_elem + pas])
    print "Calibrating rule: " + rule_name
    # Track the rule being adjusted
    rule_from_tracker.append(lu_elem)
    rule_to_tracker.append(fu_elem)
    pt_tracker.append(0)
    # Find the influence values for the specified rule
    adjust_rule_values = rules[rule_name]
    # Fix the requisite parameter values.
    old_y0 = adjust_rule_values[0]
    y1 = adjust_rule_values[1]
    y2 = adjust_rule_values[2]
    xe = adjust_rule_values[3]
    # Apply the Golden Section Search to optimise the output.
    a = min_cp
    b = max_cp
    while (abs(a - b) > gss_cp_tol):
        # Calculate the interval values.
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        # Set the neighbourhood rule to the lower interval values.
        set_lp_rule(project_file, fu_elem, lu_elem, c, y1, y2, xe)
        # Run the model, evaluate performance.
        ksim_c = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_c[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_c = sum(ksim_c) / len(ksim_c)
        # Now set the neighbourhood rule to the higher interval value.
        set_lp_rule(project_file, fu_elem, lu_elem, d, y1, y2, xe)
        # Run the model, evaluate performance.
        ksim_d = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_d[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_d = sum(ksim_d) / len(ksim_d)
        # Add step here to break out if not ksim not changing.
        if ksim_avg == ksim_c and ksim_avg == ksim_d:
            break
        else:
            # Evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            if ksim_c > ksim_d:
                b = d
            else:
                a = c
    # Set the final value, input into the model.
    if ksim_avg == ksim_c and ksim_avg == ksim_d:
        rules[rule_name][0] = old_y0
        set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
        ksim_log.append(ksim_avg)
        clu_log.append(clu_avg)
        pt_value.append(old_y0)
    else:
        final_value = (a + b) / 2
        rules[rule_name][0] = final_value
        set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
        # Evaluate the final value.
        ksim_final = [0] * max_runs
        clu_final = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_final[run_count] = ksim(omap, amap, smap, mask)
            clu_final[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                           luc, pas, act,
                                                           luc_count)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine whether to change rule or not based on result.
        ksim_avg_new = sum(ksim_final)/len(ksim_final)
        clu_avg_new = sum(clu_final)/len(clu_final)
        if ksim_avg_new < ksim_avg:
            # Conditional. If worse performance, reset.
            rules[rule_name][0] = old_y0
            set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)
            pt_value.append(old_y0)
        else:
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
            ksim_avg = ksim_avg_new
            clu_avg = clu_avg_new
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)




# 2d. Fine tuning of interaction tails. ---
print "Fine tuning the interaction tails"
# Evaluate the socio-economic groups to determine the order for
# interactive rules. First, allocate a numerical value depending on
# class converted to.
# Rank and determine order. Start by initialising an array to track order.
temp_ir_order = np.zeros(shape=(luc, act))
for i in range(0, luc):
    for j in range(0, act):
        if i == j + pas:
            # Skip if an inertia point.
            pass
        elif att_rules[i, j] == 0:
            # Skip if no interactive rule.
            pass
        else:
            # Evaluate the conversions, assign rank based on class that is
            # being converted to.
            if socio_eco_group[j + pas] == "Residential":
                temp_ir_order[i, j] = 4
            elif socio_eco_group[j + pas] == "Work":
                temp_ir_order[i, j] = 3
            elif socio_eco_group[j + pas] == "Recreation":
                temp_ir_order[i, j] = 2
            else:
                temp_ir_order[i, j] = 1
# Rank and determine order. Start by initialising an array to track order.
ir_order = np.zeros(shape=(luc, act))
# Iterate through and assign rank. Start at 1.
ranking = 1
iterator = 1
while iterator < 5:
    # Iterate through the ranks.
    for i in range(0, luc):
        for j in range(0, act):
            if temp_ir_order[i, j] == iterator:
                # If the equivalent value is found, assign current
                # ranking value.
                ir_order[i, j] = ranking
                # Add one to the ranking value.
                ranking = ranking + 1
    # Add one to the iterator.
    iterator = iterator + 1
# Apply the Golden Section Search to interactive tails.
while np.sum(ir_order) > 0:
    # Identify max value, corresponding to parameter that will be
    # tested.
    par_index = unravel_index(ir_order.argmax(), ir_order.shape)
    # Set the corresponding value in the list to 0.
    ir_order[par_index[0], par_index[1]] = 0
    # Specify the function and land-use element index values for
    # appending the geoproject file.
    lu_elem = par_index[0]
    fu_elem = par_index[1]
    # Specify the rule for tuning.
    rule_name = ("from " + luc_names[lu_elem] +
                 " to " + luc_names[fu_elem + pas])
    print "Calibrating rule: " + rule_name
    # Track the rule being adjusted
    rule_from_tracker.append(lu_elem)
    rule_to_tracker.append(fu_elem)
    pt_tracker.append(1)
    # Find the influence values for the specified rule
    adjust_rule_values = rules[rule_name]
    # Fix the requisite parameter values.
    y0 = adjust_rule_values[0]
    old_y1 = adjust_rule_values[1]
    old_y2 = adjust_rule_values[2]
    xe = adjust_rule_values[3]
    # Apply the Golden Section Search to optimise the output.
    a = min_ct
    b = max_ct
    while (abs(a - b) > gss_ct_tol):
        # Calculate the interval values.
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        # Set the neighbourhood rule to the lower interval values.
        c2 = 0.1*c
        set_lp_rule(project_file, fu_elem, lu_elem, y0, c, c2, xe)
        # Run the model, evaluate performance.
        ksim_c = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_c[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_c = sum(ksim_c) / len(ksim_c)
        # Now set the neighbourhood rule to the higher interval value.
        d2 = 0.1 * d
        set_lp_rule(project_file, fu_elem, lu_elem, y0, d, d2, xe)
        # Run the model, evaluate performance.
        ksim_d = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM
            ksim_d[run_count] = ksim(omap, amap, smap, mask)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine the average output
        ksim_d = sum(ksim_d) / len(ksim_d)
        # Add step here to break out if not ksim not changing.
        if ksim_avg == ksim_c and ksim_avg == ksim_d:
            break
        else:
            # Evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            if ksim_c > ksim_d:
                b = d
            else:
                a = c
    # Set the final value, input into the model.
    if ksim_avg == ksim_c and ksim_avg == ksim_d:
        rules[rule_name][1] = old_y1
        rules[rule_name][2] = old_y2
        set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1,
                    old_y2, xe)
        ksim_log.append(ksim_avg)
        clu_log.append(clu_avg)
        pt_value.append(old_y1)
    else:
        final_value = (a + b) / 2
        final_value_2 = final_value*0.1
        rules[rule_name][1] = final_value
        rules[rule_name][2] = final_value_2
        set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                    final_value_2, xe)
        # Evaluate the final value.
        ksim_final = [0] * max_runs
        clu_final = [0] * max_runs
        run_count = 0
        while run_count < max_runs:
            # Set the random seed
            rseed = base_seed + run_count
            set_rand(project_file, rseed)
            # Run Metronamica to generate output
            run_metro(project_file, log_file, working_directory,
                      geo_cmd)
            # Read in output
            smap = read_map(smap_path)
            # Determine KSIM & clumpiness error.
            ksim_final[run_count] = ksim(omap, amap, smap, mask)
            clu_final[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                           luc, pas, act,
                                                           luc_count)
            # Add one to iterator to prevent infinite loop!
            run_count = run_count + 1
        # Determine whether to change rule or not based on result.
        ksim_avg_new = sum(ksim_final)/len(ksim_final)
        clu_avg_new = sum(clu_final)/len(clu_final)
        if ksim_avg_new < ksim_avg:
            # Conditional. If worse performance, reset.
            rules[rule_name][1] = old_y1
            rules[rule_name][2] = old_y2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1,
                        old_y2, xe)
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)
            pt_value.append(old_y1)
        else:
            pt_value.append(final_value)
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            ksim_avg = ksim_avg_new
            clu_avg = clu_avg_new
            ksim_log.append(ksim_avg)
            clu_log.append(clu_avg)

# Save the final rules.
final_rules_file = (output_path + case_study +
                               "\\Rules\\final_rules.csv")
store = [0]*6
with open(final_rules_file, "wb") as csv_file:
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

# Save the output log.
log_output_file = (output_path + case_study + "\\fine_tuning_log.csv")
store = [0] * 6
with open(log_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    # Save and write a header line to the file
    values = ["From", "To", "Distance", "Value", "KSIM", "AWACE"]
    writer.writerow(values)
    # Now write the log
    for i in range(0, len(ksim_log)):
        store[0] = luc_names[rule_from_tracker[i]]
        store[1] = luc_names[rule_to_tracker[i] + pas]
        store[2] = pt_tracker[i]
        store[3] = pt_value[i]
        store[4] = ksim_log[i]
        store[5] = clu_log[i]
        writer.writerow(store)

# Fine tuning completed!
