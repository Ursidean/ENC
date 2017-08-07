"""
Conduct a sensitivity analysis of the ratio between points
"""

import numpy as np
from considered_distances import considered_distances
from read_map import read_map
from enrichment_factor import ef
import math
from contingency_table import contingency_table
from set_rand import set_rand
from set_NR import set_lp_rule
from kappa import kappa
from kappa import ksim
from clumpy_module import clumpiness_index
from area_weighted_clu import area_weighted_clu_error
from run_metro import run_metro
import csv

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\ENC\\"
# Specify the command line version of Geonamica.
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
# Specify the fine-tuning calibration method parameters, the base random seed, the
# maximum number of simulation runs, and the golden section search tolerance.
max_runs = 50
base_seed = 1000
# Specify the log files.
log_file_2000 = base_path + "LogSettings.xml"
log_file_2006 = base_path + "LogSettings2006.xml"


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
# Specify the fixed transformation parameters.

theta_i1 = 0.037
theta_i2 = theta_i1*0.1

theta_cp = 0.010

theta_c1 = 0.017
theta_c2 = theta_c1*0.1

# Calculate all relevant parameter values.
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

# Initialise a dictionary to store the metrics calculated.
metrics = {}

# Specify a list of case studies to analyse.
case_studies = ["Rome"]#["Berlin", "Budapest", "Madrid", "Rome"]

# Iterate through the case studies for testing.
for a in range(0, len(case_studies)):
    # Specify the selected case study.
    case_study = case_studies[a]
    # Terminal feedback to user.
    print "Testing case study: " + case_study
    # Set the paths to the directories and relevant data
    data_path = base_path + "Example_case_study_data\\"
    output_path = base_path + "Example_case_study_output\\" 
    map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
    map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
    map3_path = data_path + case_study + "\\" + case_study.lower() + "_2006.asc"
    mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
    # Specify the working directory, the folder containing the Metronamica
    # project file.
    working_directory_2000 = ("C:\\Users\\charl\\OneDrive\\Documents"
                              "\\Geonamica\\Metronamica\\"
                              + case_study + "\\")
    working_directory_2006 = ("C:\\Users\\charl\\OneDrive\\Documents"
                              "\\Geonamica\\Metronamica\\"
                              + case_study + "_2006\\")
    # Specify the project file names.
    project_file_2000 = working_directory_2000 + case_study + ".geoproj"
    project_file_2006 = working_directory_2006 + case_study + "_2006.geoproj"

    # Specify the paths to the simulated output maps.
    smap_path_2000 = ("C:\\Users\\charl\\OneDrive\\Documents"
                      "\\Geonamica\\Metronamica\\"
                      + case_study + "\\Log\\Land_use\\"
                      "Land use map_2000-Jan-01 00_00_00.rst")
    smap_path_2006 = ("C:\\Users\\charl\\OneDrive\\Documents"
                      "\\Geonamica\\Metronamica\\"
                      + case_study + "_2006\\Log\\Land_use\\"
                      "Land use map_2006-Jan-01 00_00_00.rst")

    # Read in the data maps and mask.
    map1 = read_map(map1_path)
    map2 = read_map(map2_path)
    map3 = read_map(map3_path)
    mask = read_map(mask_path)
    # Analyse the input maps for evaluation purposes
    map_dimensions = np.shape(map1)
    rows = map_dimensions[0]
    cols = map_dimensions[1]
    # Generate the Enrichment Factor and contingency table.
    data_ef = ef(luc, max_distance, cdl, cd, N, map1, map2, mask, rows, cols)
    # Log scale the enrichment factor values.
    log_data_ef = np.zeros(shape=(max_distance, luc, luc))
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_distance):
                if data_ef[c, p, q] == 0:
                    log_data_ef[c, p, q] = -9999
                else:
                    log_data_ef[c, p, q] = math.log(data_ef[c, p, q], 10)
    # Generate the contingency table using the module, 'contingency_table'
    cont_table = contingency_table(map1, map2, mask, luc, rows, cols)
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
            set_lp_rule(project_file_2000, fu_elem, lu_elem, y0, y1, y2, xe)
            set_lp_rule(project_file_2006, fu_elem, lu_elem, y0, y1, y2, xe)
    run_count = 0
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
        key1 = case_study + "-kappa-2000-it-" + str(seed)
        metrics[key1] = it_kappa
        key2 = case_study + "-ksim-2000-it-" + str(seed)
        metrics[key2] = it_ksim
        key3 = case_study + "-clu-2000-it-" + str(seed)
        metrics[key3] = it_clu
        # Perform analysis for validation period (2000-2006).

        set_rand(project_file_2006, seed)
        run_metro(project_file_2006, log_file_2006, working_directory_2006,
                    geo_cmd)
        smap = read_map(smap_path_2006)
        it_kappa = kappa(map3, smap, mask)
        it_ksim = ksim(map2, map3, smap, mask)
        it_clu = clumpiness_index(smap, mask, luc)
        key1 = case_study + "-kappa-2006-it-" + str(seed)
        metrics[key1] = it_kappa
        key2 = case_study + "-ksim-2006-it-" + str(seed)
        metrics[key2] = it_ksim
        key3 = case_study + "-clu-2006-it-" + str(seed)
        metrics[key3] = it_clu

        # Add 1 to iterator
        run_count = run_count + 1

# Write the output to a .csv file.
metrics_output_file = output_path + case_study + "\\meta_metrics.csv"
store = [0]*(5 + luc)
with open (metrics_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["Case study", "Year", "Seed", "Kappa", "KSIM"] + luc_names
    writer.writerow(values)
    for a in range(0, len(case_studies)):
        case_study = case_studies[a]
        for x in range(0, max_runs):
            seed = base_seed + x
            store[0] = case_study
            store[1] = "2000"

            store[2] = str(x)
            key1 = case_study + "-kappa-2000-it-" + str(seed)
            key2 = case_study + "-ksim-2000-it-" + str(seed)
            key3 = case_study + "-clu-2000-it-" + str(seed)
            store[3] = metrics[key1]
            store[4] = metrics[key2]
            for c in range(0, luc):
                store[5 + c] = metrics[key3][c]
            writer.writerow(store)

        for x in range(0, max_runs):
            seed = base_seed + x
            store[0] = case_study
            store[1] = "2006"
            store[2] = str(x)
            key1 = case_study + "-kappa-2006-it-" + str(seed)
            key2 = case_study + "-ksim-2006-it-" + str(seed)
            key3 = case_study + "-clu-2006-it-" + str(seed)
            store[3] = metrics[key1]
            store[4] = metrics[key2]
            for c in range(0, luc):
                store[5 + c] = metrics[key3][c]
            writer.writerow(store)

                
import winsound
freq=2500
dur=1000
winsound.Beep(freq,dur)
