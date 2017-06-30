"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine tuning of the neighbourhood rule ratios to improve the agreement with the
observed enrichment factor values.
"""

import time
import math
import numpy as np
from read_map import read_map
from considered_distances import considered_distances
from contingency_table import contingency_table
from enrichment_factor import ef
from log_scale_ef import log_scale_ef
from set_NR import set_lp_rule
from multi_sims_ef import multi_sims_ef
from kappa import ksim
import csv
from ef_ssr_calc import ef_ssr

start_time = time.time()

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Berlin"
# Set the paths to the relevant directories
data_path = base_path + "Example_case_study_data\\"
output_path = base_path + "Example_case_study_output\\"
working_directory = ("C:\\Users\\charl\\OneDrive\\Documents"
                     "\\Geonamica\\Metronamica\\" +
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
# Specify the fine-tuning calibration method parameters, the base random seed,
# the maximum number of simulation runs, the golden section search tolerance,
# and an error value for when the sum of square residuals can't be evaluated.
# (default value of 5.0).
max_runs = 3
base_seed = 1000
no_ssr_eval_error = 5.0
log_base = 10
# Specify the golden ratio (approx.).
gr = (math.sqrt(5) + 1)/2

# Band settings.
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

# Specify the Inertia values.
high_inertia = 1000.0
mid_inertia = 500.0
low_inertia = 250.0
# Specify the initial transformation parameters.
theta_cp = 0.025
theta_i1 = 0.05
theta_i2 = theta_i1*0.1
theta_c1 = 0.005
theta_c2 = theta_c1*0.1

# Calculate all relevant neighbourhood rule parameter values.
# Set the conversion parameter values.
high_conversion = high_inertia * theta_cp
mid_conversion = mid_inertia * theta_cp
low_conversion = low_inertia * theta_cp
# Set the distance 1 tail values for self-influence rules.
d1_high_si_value = high_inertia * theta_i1
d1_mid_si_value = mid_inertia * theta_i1
d1_low_si_value = low_inertia * theta_i1
# Set the distance 2 tail values for self-influence rules.
d2_high_si_value = high_inertia * theta_i2
d2_mid_si_value = mid_inertia * theta_i2
d2_low_si_value = low_inertia * theta_i2
# Set the distance 1 tail values for interaction rules.
d1_high_co_value = high_inertia * theta_c1
d1_mid_co_value = mid_inertia * theta_c1
d1_low_co_value = low_inertia * theta_c1
# Set the distance 2 tail value for interaction rules.
d2_high_co_value = high_inertia * theta_c2
d2_mid_co_value = mid_inertia * theta_c2
d2_low_co_value = low_inertia * theta_c2

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

# Golden Section Search (GSS) variable initialisation.
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
# Generate the initial Sum of Square Residuals (SSR) values
dummy = ef_ssr(log_data_ef, log_sims_ef, luc, act, no_ssr_eval_error)
ssr_matrix = dummy[0]
ssr_sum = dummy[1]
ssr_sum_log.append(ssr_sum)
smap = read_map(smap_path)
initial_ksim = ksim(omap, amap, smap, mask)
ksim_log.append(initial_ksim)

# Begin the fine-tuning of the ratio between different points.
# First, fine-tune the inertia tail ratio.
theta_i1_a = 0.0
theta_i1_b = 0.5
tol_it = 0.05
while abs(theta_i1_b - theta_i1_a) > tol_it:
    theta_i1_c = theta_i1_b - (theta_i1_b - theta_i1_a)/gr
    theta_i1_d = theta_i1_a + (theta_i1_b - theta_i1_a)/gr
    theta_i2_c = theta_i1_c*0.1
    theta_i2_d = theta_i1_d * 0.1

    # Re-calculate the neighbourhood rule values for the lower bound.
    # Set the distance 1 tail values for self-influence rules.
    d1_high_si_value = high_inertia * theta_i1_c
    d1_mid_si_value = mid_inertia * theta_i1_c
    d1_low_si_value = low_inertia * theta_i1_c
    # Set the distance 2 tail values for self-influence rules.
    d2_high_si_value = high_inertia * theta_i2_c
    d2_mid_si_value = mid_inertia * theta_i2_c
    d2_low_si_value = low_inertia * theta_i2_c
    # Update the rule set
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
    # Calculate the simulated enrichment factors, log scale.
    c_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path,project_file,
                              log_file,working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_c_sims_ef = log_scale_ef(c_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_c = ef_ssr(log_data_ef, log_c_sims_ef, luc, act, no_ssr_eval_error)
    f_c = ssr_c[1]

    # Re-calculate the neighbourhood rule values for the higher bound.
    # Set the distance 1 tail values for self-influence rules.
    d1_high_si_value = high_inertia * theta_i1_d
    d1_mid_si_value = mid_inertia * theta_i1_d
    d1_low_si_value = low_inertia * theta_i1_d
    # Set the distance 2 tail values for self-influence rules.
    d2_high_si_value = high_inertia * theta_i2_d
    d2_mid_si_value = mid_inertia * theta_i2_d
    d2_low_si_value = low_inertia * theta_i2_d
    # Update the rule set
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
    # Calculate the simulated enrichment factors, log scale.
    d_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                              log_file, working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_d_sims_ef = log_scale_ef(d_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_d = ef_ssr(log_data_ef, log_d_sims_ef, luc, act, no_ssr_eval_error)
    f_d = ssr_d[1]
    # Evaluate the result
    if f_c < f_d:
        theta_i1_b = theta_i1_d
    else:
        theta_i1_a = theta_i1_c
# Process the final evaluation, and update the rule.
if f_c <= f_d:
    theta_i1 = theta_i1_c
else:
    theta_i1 = theta_i1_d

theta_i2 = theta_i1*0.1

# Print the inertia tail ratio.
print theta_i1

# Set the distance 1 tail values for self-influence rules.
d1_high_si_value = high_inertia * theta_i1
d1_mid_si_value = mid_inertia * theta_i1
d1_low_si_value = low_inertia * theta_i1
# Set the distance 2 tail values for self-influence rules.
d2_high_si_value = high_inertia * theta_i2
d2_mid_si_value = mid_inertia * theta_i2
d2_low_si_value = low_inertia * theta_i2
# Update the rule set
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

# Now tune the conversion points.
theta_cp_a = 0.0
theta_cp_b = 0.2
tol_cp = 0.02
while abs(theta_cp_a - theta_cp_b) > tol_it:
    theta_cp_c = theta_cp_b - (theta_cp_b - theta_cp_a)/gr
    theta_cp_d = theta_cp_a + (theta_cp_b - theta_cp_a)/gr

    # Re-calculate the neighbourhood rule values for the lower bound.
    high_conversion = high_inertia * theta_cp_c
    mid_conversion = mid_inertia * theta_cp_c
    low_conversion = low_inertia * theta_cp_c
    # Update the rule set.
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
    # Calculate the simulated enrichment factors, log scale.
    c_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                              log_file, working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_c_sims_ef = log_scale_ef(c_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_c = ef_ssr(log_data_ef, log_c_sims_ef, luc, act, no_ssr_eval_error)
    f_c = ssr_c[1]

    # Re-calculate the neighbourhood rule values for the higher bound.
    high_conversion = high_inertia * theta_cp_d
    mid_conversion = mid_inertia * theta_cp_d
    low_conversion = low_inertia * theta_cp_d
    # Update the rule set.
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
    # Calculate the simulated enrichment factors, log scale.
    d_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                              log_file, working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_d_sims_ef = log_scale_ef(d_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_d = ef_ssr(log_data_ef, log_d_sims_ef, luc, act, no_ssr_eval_error)
    f_d = ssr_d[1]
    # Evaluate the result
    if f_c < f_d:
        theta_cp_b = theta_cp_d
    else:
        theta_cp_a = theta_cp_c
# Process the final evaluation, and update the rule.
if f_c <= f_d:
    theta_cp = theta_cp_c
else:
    theta_cp = theta_cp_d

# Print the conversion point ratio.
print theta_cp

# Set the conversion parameter values.
high_conversion = high_inertia * theta_cp
mid_conversion = mid_inertia * theta_cp
low_conversion = low_inertia * theta_cp
# Update the rule set.
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

# Finally, fine-tune the conversion tail ratio.
theta_c1_a = 0.0
theta_c1_b = 0.1
tol_ct = 0.01
while abs(theta_c1_b - theta_c1_a) > tol_ct:
    theta_c1_c = theta_c1_b - (theta_c1_b - theta_c1_a)/gr
    theta_c1_d = theta_c1_a + (theta_c1_b - theta_c1_a)/gr
    theta_c2_c = theta_c1_c*0.1
    theta_c2_d = theta_c1_d * 0.1

    # Set the distance 1 tail values for interaction rules.
    d1_high_co_value = high_inertia * theta_c1_c
    d1_mid_co_value = mid_inertia * theta_c1_c
    d1_low_co_value = low_inertia * theta_c1_c
    # Set the distance 2 tail value for interaction rules.
    d2_high_co_value = high_inertia * theta_c2_c
    d2_mid_co_value = mid_inertia * theta_c2_c
    d2_low_co_value = low_inertia * theta_c2_c
    # Update the rule set
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, pass.
            if i + pas == j:
                pass
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
    # Calculate the simulated enrichment factors, log scale.
    c_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path,project_file,
                              log_file,working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_c_sims_ef = log_scale_ef(c_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_c = ef_ssr(log_data_ef, log_c_sims_ef, luc, act, no_ssr_eval_error)
    f_c = ssr_c[1]

    # Re-calculate the neighbourhood rule values for the higher bound.
    # Set the distance 1 tail values for interaction rules.
    d1_high_co_value = high_inertia * theta_c1_d
    d1_mid_co_value = mid_inertia * theta_c1_d
    d1_low_co_value = low_inertia * theta_c1_d
    # Set the distance 2 tail value for interaction rules.
    d2_high_co_value = high_inertia * theta_c2_d
    d2_mid_co_value = mid_inertia * theta_c2_d
    d2_low_co_value = low_inertia * theta_c2_d
    # Update the rule set
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, pass.
            if i + pas == j:
                pass
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
    # Calculate the simulated enrichment factors, log scale.
    d_sims_ef = multi_sims_ef(max_runs, base_seed, smap_path, project_file,
                              log_file, working_directory, geo_cmd, luc,
                              max_distance, cdl, cd, N, omap, mask,
                              rows, cols)
    log_d_sims_ef = log_scale_ef(d_sims_ef, log_base, luc, act, pas,
                                 max_distance)
    # Calculate the sum of square residuals for the lower interval.
    ssr_d = ef_ssr(log_data_ef, log_d_sims_ef, luc, act, no_ssr_eval_error)
    f_d = ssr_d[1]
    # Evaluate the result
    if f_c < f_d:
        theta_c1_b = theta_c1_d
    else:
        theta_c1_a = theta_c1_c
# Process the final evaluation, and update the rule.
if f_c <= f_d:
    theta_c1 = theta_c1_c
else:
    theta_c1 = theta_c1_d
theta_c2 = theta_c1*0.1

# Print the conversion tail ratio.
print theta_c1

# Set the distance 1 tail values for interaction rules.
d1_high_co_value = high_inertia * theta_c1
d1_mid_co_value = mid_inertia * theta_c1
d1_low_co_value = low_inertia * theta_c1
# Set the distance 2 tail value for interaction rules.
d2_high_co_value = high_inertia * theta_c2
d2_mid_co_value = mid_inertia * theta_c2
d2_low_co_value = low_inertia * theta_c2
# Update the rule set
for i in range(0, act):
    for j in range(0, luc):
        # Specify the neighbourhood rule key.
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        # If a self-influence rule, pass.
        if i + pas == j:
            pass
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

# Save the final set of rules.
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
# Determine the time taken.
end_time = time.time()
total_time = end_time - start_time
print "Time taken:" + str(total_time)
# Calibration completed!
