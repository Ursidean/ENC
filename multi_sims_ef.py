"""
This module averages enrichment factor values calculated for multiple model simulations
"""

import numpy as np
from enrichment_factor import ef
from run_metro import run_metro
from set_rand import set_rand
from read_map import read_map


def multi_sims_ef(max_runs, base_seed, smap_path, project_file, log_file,
                  working_directory, geo_cmd, luc, max_d, cdl, cd, N,
                  omap, mask, rows, cols):
    # Initialise an run counter to track the number of iterations performed.
    run_count = 0
    # Initialise a dictionary to store simulated enrichment factor values.
    sim_ef_store = {}
    while run_count < max_runs:
        # Set the dictionary key as the run count.
        key = run_count
        # Generate the random seed value
        random_seed = base_seed + run_count
        # Set the random seed for the metronamica geoproject file
        set_rand(project_file, random_seed)
        # Run the model to generate a simulated output, read in the output map.
        run_metro(project_file, log_file, working_directory, geo_cmd)
        smap = read_map(smap_path)
        # Calculate the enrichment factor for the simulation
        sim_ef = ef(luc, max_d, cdl, cd, N, omap, smap, mask, rows, cols)
        # Store in the dictionary, add 1 to the iterator.
        sim_ef_store[key] = sim_ef
        run_count = run_count + 1
    # Calculate the average enrichment factor across the random seeds
    sim_ef_average = np.zeros(shape=(max_d, luc, luc))
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_d):
                # Initialise a dummy storage array to store values for
                # evaluation.
                store = [0]*max_runs
                for i in range(0, max_runs):
                    store[i] = sim_ef_store[i][c, p, q]
                sim_ef_average[c, p, q] = (sum(store))/(len(store))
    return sim_ef_average
