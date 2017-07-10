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
from set_rand import set_rand
from run_metro import run_metro
from kappa import ksim
from AACE import AACE
import csv

start_time = time.time()

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Select an example case study application. Specify the name below:
case_study = "Rome"
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
# Determine the size of each land-use class.
luc_size = [0]*luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_size[amap[i, j]] = luc_size[amap[i, j]] + 1

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
# Specify the calibration method parameters.
max_runs = 10
base_seed = 1000

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
# Initialise a dictionary to store metric values.
metrics = {}
# Input the rules to be analysed. Start by initialising a dictionary for
# storage.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, 5]
# Set the step-size and bounds for each parameter
low_bound = 0
high_bound = 21
theta_it_step = 0.01
theta_cp_step = 0.01
theta_ct_step = 0.001
# First, enumerate across the range of ratio values for inertia tail.
print "Enumerating self-influence tail values"
# Set the fixed components.
theta_cp = 0.025
theta_c1 = 0.005
theta_c2 = theta_c1 * 0.1
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
# Input into model.






for b in range(low_bound, high_bound, 1):
    # Set the variable parameter value. This will change depending on
    # what is tested. Hard coded for convenience.
    varied_parameter = "theta_it"
    theta_i1 = float(b) * theta_it_step
    theta_i2 = 0.1 * theta_i1
    print theta_i1
    # Set the distance 1 tail values for self-influence rules.
    d1_high_si_value = high_inertia * theta_i1
    d1_mid_si_value = mid_inertia * theta_i1
    d1_low_si_value = low_inertia * theta_i1
    # Set the distance 2 tail values for self-influence rules.
    d2_high_si_value = high_inertia * theta_i2
    d2_mid_si_value = mid_inertia * theta_i2
    d2_low_si_value = low_inertia * theta_i2
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
    # Run the simulation model, calculate the metric values.
    run_count = 0
    # Now generate the simulated output and track the results.
    while run_count < max_runs:
        # Perform analysis for the calibration period (1990-2000).
        seed = base_seed + run_count
        set_rand(project_file, seed)
        run_metro(project_file, log_file, working_directory,
                  geo_cmd)
        smap = read_map(smap_path)
        it_ksim = ksim(omap, amap, smap, mask)
        it_clu = AACE(amap, smap, mask, pas, act, luc, luc_size)
        key1 = (
            case_study + "-" + varied_parameter + "-" + str(b)
            + "-ksim-2000-it-" + str(seed)
        )
        metrics[key1] = it_ksim
        key2 = (
            case_study + "-" + varied_parameter + "-" + str(b)
            + "-clu-2000-it-" + str(seed)
        )
        metrics[key2] = it_clu
        # Add one to iterator
        run_count = run_count + 1
    # Select a value?
    # Code?
# Write the output to a .csv file.
metrics_output_file = (
    output_path + case_study + "\\meta_parameter_theta_it.csv"
)
# Store the output metric values.
store = [0]*4
with open (metrics_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Var Par Value", "Seed", "KSIM", "AWACE"]
    writer.writerow(values)
    # First, write the values for inertia tail values.
    for b in range(low_bound, high_bound, 1):
        varied_parameter = "theta_it"
        theta_i1 = float(b) * theta_it_step
        store[0] = str(theta_i1)
        for x in range(0, max_runs):
            seed = base_seed + x
            store[1] = str(x)
            key1 = (
                case_study + "-" + varied_parameter + "-" + str(b)
                + "-ksim-2000-it-" + str(seed)
            )
            key2 = (
                case_study + "-" + varied_parameter + "-" + str(b)
                + "-clu-2000-it-" + str(seed)
            )
            store[2] = metrics[key1]
            store[3] = metrics[key2]
            writer.writerow(store)

# Determine the time taken.
end_time = time.time()
total_time = end_time - start_time
print "Time taken:" + str(total_time)
# Enumeration completed!
