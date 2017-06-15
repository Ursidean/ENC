"""
This module calculates the Sum of Square Residuals for 2 log-scaled matrices
of enrichment factor values for different distances(simulation and data).
"""

import numpy as np


def log_based_ssr_calc(data_ef_line1, sim_ef_line2):
    # Calculate the error bewteen the two enrichment factor values
    length = len(data_ef_line1)
    sq_errors = [0]*length
    for i in range(0, length):
        # If both return a non-evlauation code, skip calculation
        if data_ef_line1[i] < -1000 and sim_ef_line2[i] < -1000:
            pass
        # If one is evaluated and one is not, return an error code -9999.
        # Alternatively could return a high error value
        elif data_ef_line1[i] < -1000 and sim_ef_line2[i] > -1000:
            return -9999
        elif data_ef_line1[i] > -1000 and sim_ef_line2[i] < -1000:
            return -9999
        else:
            sq_errors[i] = (data_ef_line1[i] - sim_ef_line2[i])**2
    ssr = sum(sq_errors)
    return ssr


def ef_ssr(log_ef_1, log_ef_2, luc, act, no_ssr_eval_error):
    ssr_values = np.zeros(shape=(luc, act))
    ssr_sum = 0

    for i in range(0, luc):
        for j in range(0, act):
            ef_line1 = log_ef_1[:, i, j]
            ef_line2 = log_ef_2[:, i, j]
            ef_line1.tolist()
            ef_line2.tolist()
            x = log_based_ssr_calc(ef_line1, ef_line2)
            if x == -9999:
                ssr_values[i, j] = no_ssr_eval_error
                ssr_sum = ssr_sum + no_ssr_eval_error
            else:
                ssr_values[i, j] = x
                ssr_sum = ssr_sum + x
    return ssr_values, ssr_sum
