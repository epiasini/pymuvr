"""
This is a standalone benchmark script that can be used to compare
execution times between pymuvr and spykeutils.
"""
import numpy as np
import time

import pymuvr

import quantities as pq
import spykeutils.spike_train_generation as stg
import spykeutils.spike_train_metrics as stm
    
def generate_spike_trains(n_observations, n_cells, rate, tstop):
    sutils_units = {}
    pymuvr_observations = []
    for unit in range(n_cells):
        sutils_units[unit] = []
        for ob in range(n_observations):
            sutils_units[unit].append(stg.gen_homogeneous_poisson(rate * pq.Hz, t_stop=tstop * pq.s))
    # observation 1 is identical to observation 0 for all the cells.
    for unit in range(n_cells):
        sutils_units[unit][1] = sutils_units[unit][0]

    for ob in range(n_observations):
        pymuvr_observations.append([])
        for unit in range(n_cells):
            pymuvr_observations[ob].append(sutils_units[unit][ob].tolist())

    return sutils_units, pymuvr_observations

def main():
    n_observations = 5
    n_cells = 50
    rate = 30
    tstop = 3
    cos = 0.1
    tau = 0.012

    sutils_units, pymuvr_observations = generate_spike_trains(n_observations, n_cells, rate, tstop)
    sutils_start = time.clock()
    sutils_d = stm.van_rossum_multiunit_dist(sutils_units, weighting=cos, tau=tau)
    sutils_stop = time.clock()
    pymuvr_start = time.clock()
    pymuvr_d = pymuvr.square_distance_matrix(pymuvr_observations, cos, tau)
    pymuvr_stop = time.clock()

    print("\nspykeutils distance matrix:\n{0}\n\n".format(sutils_d))
    print("\npymuvr distance matrix:\n{0}\n\n".format(pymuvr_d))
    print("Time elapsed: pymuvr {0}s, spykeutils {1}s".format(pymuvr_stop-pymuvr_start, sutils_stop-sutils_start))


if __name__ == "__main__":
    main()
