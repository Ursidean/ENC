"""
This component of the Empirical Neighbourhood Calibration method sets the initial
neighbourhood rule values
"""
import gdal
import numpy as np
import csv
import math
from considered_distances import considered_distances
from enrichment_factor import ef
from contingency_table import contingency_table


# Specify the initialisation settings.
# Specify the enrichment factor band levels, for setting the tail values.
high_ef = 1.0
mid_ef = 0.5
# Specify the bands for the different levels of inertia.
high_inertia_band = 0.95
mid_inertia_band = 0.90
# Specify the different inertia values taken based on band.
high_inertia = 1000
mid_inertia = 800
low_inertia = 500
# Specify the tail value for self-influence rules at a distance of 1.
d1_high_si_value = 100
d1_mid_si_value = 50
d1_low_si_value = 20
# Specify the tail value for self-influence rules at a distance of 2.
d2_high_si_value = 10
d2_mid_si_value = 5
d2_low_si_value = 2
# Specify the bands for the different levels of conversion.
high_conversion_band = 0.5
mid_conversion_band = 0.1
low_conversion_band = 0.02
# Specify the different conversion values taken based on band
high_conversion = 50
mid_conversion = 20
low_conversion = 5
# Specify the tail value for interaction rules at a distance of 1.
d1_high_co_value = 10
d1_mid_co_value = 5
d1_low_co_value = 2
# Specify the tail value for interation rules at a distance of 2.
d2_high_co_value = 1
d2_mid_co_value = 0.5
d2_low_co_value = 0.2

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
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

# Calculate the enrichment factor values for the data for initialisation setting,
# using the module, 'ef'.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
# Log scale the enrichmet factor values
log_data_ef = np.zeros(shape =(max_distance, luc, luc))
for p in range(0, luc):
    for q in range(0, luc):
        for c in range(0, max_distance):
            if data_ef[c, p, q] == 0:
                log_data_ef[c, p, q] = -9999
            else:
                log_data_ef[c, p, q] = math.log(data_ef[c, p, q], 10)
# Generate the contingency table using the module, 'contingency_table'
cont_table = contingency_table(omap, amap, mask, luc, rows, cols)

# Generate the initial rule set by evaluating the contingency table and enrichment
# factor values. Store the rules in a dictionary, indexed as from ... to ...
# The rules are formatted as: Influence values at distance 0, 1, 2, and the
# distance when influence is 0.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, max_distance]
# Load the attraction rules file
att_rule_file = output_path + case_study + "\\Rules\\att_rules.txt"
att_rules = np.loadtxt(att_rule_file)
# Evaluate the rates of inertia and conversion
ic_rates = np.zeros(shape=(luc, luc))
for i in range(0, luc):
    for j in range(0, luc):
        if i == j:
            if cont_table[i, luc] > 0:
                ic_rates[i, j] = cont_table[i, j]/cont_table[i, luc]
        else:
            conversions = abs(float(cont_table[j, j]) - float(cont_table[luc, j]))
            if conversions > 0:
                ic_rates[i, j] = float(cont_table[i, j]) / float(conversions)
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
                        if log_data_ef[c, i + pas, j] > high_ef:
                            rules[key][c] = d1_high_si_value
                        elif log_data_ef[c, i + pas, j] > mid_ef:
                            rules[key][c] = d1_mid_si_value
                        else:
                            rules[key][c] = d1_low_si_value
                elif c == 2:
                    if att_rules[j, i] == 1:
                        if log_data_ef[c, i + pas, j] > high_ef:
                            rules[key][c] = d2_high_si_value
                        elif log_data_ef[c, i + pas, j] > mid_ef:
                            rules[key][c] = d2_mid_si_value
                        else:
                            rules[key][c] = d2_low_si_value
        # If a conversion rule, set the interactive attraction values.
        else:
            if (att_rules[j, i] == 1 and log_data_ef[1, i + pas, j] > 0 and
               log_data_ef[2, i + pas, j] > 0):
                for c in range(1, 3):
                    if c == 1:
                        if log_data_ef[c, i + pas, j] > high_ef:
                            rules[key][c] = d1_high_co_value
                        elif log_data_ef[c, i + pas, j] > mid_ef:
                            rules[key][c] = d1_mid_co_value
                        elif log_data_ef[c, i + pas, j] > 0:
                            rules[key][c] = d1_low_co_value
                    elif c == 2:
                        if log_data_ef[c, i + pas, j] > high_ef:
                            rules[key][c] = d2_high_co_value
                        elif log_data_ef[c, i + pas, j] > mid_ef:
                            rules[key][c] = d2_mid_co_value
                        elif log_data_ef[c, i + pas, j] > 0:
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
                if att_rules[j, i] == 1 and log_data_ef[c, i + pas, j] > 0:
                    rules[key][3] = c + 1

# Save the initial rules generated.
initial_rules_file = output_path + case_study + "\\Rules\\initial_rules.csv"
store = [0] * 6
with open(initial_rules_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["from", "to", "I@x=0", "I@x=1", "I@x=2", "xwI=0"]
    writer.writerow(values)
    for i in range(0, luc):
        for j in range(0, act):
            # fu_elem=i
            # lu_elem=j
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            store[0] = luc_names[i]
            store[1] = luc_names[j + pas]
            store[2] = rules[key][0]
            store[3] = rules[key][1]
            store[4] = rules[key][2]
            store[5] = rules[key][3]
            writer.writerow(store)
