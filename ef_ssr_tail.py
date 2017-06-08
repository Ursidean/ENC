"""
This module calculates the Sum of Square Residuals for 2 log-scaled matrices
of enrichment factor values for different distances(simulation and data).
"""

import numpy as np


def log_based_ssr_calc(ef_line1, ef_line2):
    length = len(ef_line1)




def ef_ssr(log_ef_1, log_ef_2, att_rules, luc, act):
    ssr_values = np.zeros(shape=(luc, act))
    ssr_sum = 0
    for i in range(0, luc):
        for j in range(0, act):
            if att_rules[i, j] == 1:
                ef_line1 = log_ef_1[:, i, j]
                ef_line2 = log_ef_2[:, i, j]
                ef_line1.tolist()
                ef_line2.tolist()
