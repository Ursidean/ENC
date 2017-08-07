"""
Evaluation and plotting of a Pareto front extracted from a set of points.
"""

import numpy as np
import matplotlib.pyplot as plt

# Generate a set of data.
a = np.array([[0.001429515, 0.001839971],
             [0.001450598, 0.001851585],
             [0.001798744, 0.001864911],
             [0.001762742, 0.001836292],
             [0.002265497, 0.001847795],
             [0.003251971, 0.001852451],
             [0.003363188, 0.001848725],
             [0.004966808, 0.001850028],
             [0.006495707, 0.00187893],
             [0.006953137, 0.001961151],
             [0.008377043, 0.002073655],
             [0.009884883, 0.002211371],
             [0.011254595, 0.002343509],
             [0.012331564, 0.002512046],
             [0.01453014, 0.002736347],
             [0.016921865, 0.00296579],
             [0.017917513, 0.003203624],
             [0.019900123, 0.003439426],
             [0.02207523, 0.003710807],
             [0.022968885, 0.003971036],
             [0.023000568, 0.004220012],
             [0.023782214, 0.004484575],
             [0.024148762, 0.004777428],
             [0.023763159, 0.005089897],
             [0.023281965, 0.005441295],
             [0.022728719, 0.005787028],
             [0.021758259, 0.006163707],
             [0.021127195, 0.006506409],
             [0.020676256, 0.006909024],
             [0.020410868, 0.007283659],
             [0.019389774, 0.007681136],
             [0.018926105, 0.008102495],
             [0.018379404, 0.008566754],
             [0.017392503, 0.009091398],
             [0.01710634, 0.009573159],
             [0.016730453, 0.010082093]])

# Pareto front extraction
# Assumed that want to maximise column 1 and minimise column 2.
# Determine the number of points to evaluate
points = a.shape[0]
# Initialise a list of True/False, where True is on the Pareto front.
is_pareto = np.ones(a.shape[0], dtype = bool)
for i in range(0, points):
    x1 = a[i, 0]
    y1 = a[i, 1]
    for j in range(0, points):
        if i == j:
            pass
        else:
            x2 = a[j, 0]
            y2 = a[j, 1]
            if x2 >= x1 and y2 <= y1:
                is_pareto[i] = 0
                break
# Extract an array of the Pareto front.
pareto_points = np.sum(is_pareto)
Pareto_front = np.zeros(shape = (pareto_points, 2))
track = 0
for i in range(0, points):
    if is_pareto[i] == 1:
        Pareto_front[track, 0] = a[i, 0]
        Pareto_front[track, 1] = a[i, 1]
        track = track + 1

# Find the mid-point solution
mid_point_index = int((float(Pareto_front.shape[0]))/2)
mid_point_solution = [Pareto_front[mid_point_index, 0], Pareto_front[mid_point_index, 1]]

# Plot the points, highlighting the Pareto front.
# Extract the pareto front points for x & y
all_x = a[:, 0]
all_y = a[:, 1]
pareto_x = Pareto_front[:, 0]
pareto_y = Pareto_front[:, 1]
plt.scatter(all_x, all_y, c="b", label = "all")
plt.scatter(pareto_x, pareto_y, c="r", label = "Pareto front")
plt.show()
