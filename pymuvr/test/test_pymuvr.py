import unittest
import random

import pymuvr

try:
    import numpy as np
    NUMPY_IS_AVAILABLE = True

    def spiketrain_to_list(spiketrain):
        return [np.double(t) for t in spiketrain]
except ImportError:
    NUMPY_IS_AVAILABLE = False

try:
    import quantities as pq
    import spykeutils.spike_train_generation as stg
    import spykeutils.spike_train_metrics as stm
    SPYKEUTILS_IS_AVAILABLE = True
except ImportError:
    SPYKEUTILS_IS_AVAILABLE = False


def simple_train(mean_isi, max_duration):
    train = []
    last_spike = 0
    while last_spike < max_duration:
        delta = random.uniform(0, 2*mean_isi)
        last_spike += delta
        train.append(last_spike)
    return train

class TestDistanceMatrix(unittest.TestCase):
    def setUp(self):
        n_observations = 10
        n_cells = 100
        mean_isi = 0.03
        max_duration = 2
        self.tau = 0.012
        self.cos = 0.5
        self.observations = [[simple_train(mean_isi, max_duration) for c in range(n_cells)] for o in range(n_observations)]
        # observation 1 is identical to observation 0 for all the cells.
        self.observations[0] = self.observations[1][:]
    def test_square_distance_matrix(self):
        d = pymuvr.square_distance_matrix(self.observations, self.cos, self.tau)
        self.assertEqual(d.shape, (len(self.observations), len(self.observations)))
    def test_distance_matrix(self):
        d = pymuvr.distance_matrix(self.observations[:3],
                                   self.observations[3:],
                                   self.cos,
                                   self.tau)
        self.assertEqual(d.shape, (3, len(self.observations)-3))
    @unittest.skipIf(not NUMPY_IS_AVAILABLE, "can't import numpy")
    def test_compare_square_and_rectangular(self):
        d_rectangular = pymuvr.distance_matrix(self.observations,
                                               self.observations,
                                               self.cos,
                                               self.tau)
        d_square = pymuvr.square_distance_matrix(self.observations,
                                                 self.cos,
                                                 self.tau)

        np.testing.assert_array_almost_equal(d_rectangular, d_square)

@unittest.skipIf(not SPYKEUTILS_IS_AVAILABLE or not NUMPY_IS_AVAILABLE,
                 "can't import spykeutils or numpy, or both")
class TestCompareWithSpykeutils(unittest.TestCase):
    def setUp(self):
        self.n_observations = 10
        self.n_cells = 20
        self.rate = 30
        self.tstop = 2
        self.cos = 0.1
        self.tau = 0.012
        self.sutils_units = {}
        self.pymuvr_observations = []
        for unit in range(self.n_cells):
            self.sutils_units[unit] = []
            for ob in range(self.n_observations):
                self.sutils_units[unit].append(stg.gen_homogeneous_poisson(self.rate * pq.Hz, t_stop=self.tstop * pq.s))
        # observation 1 is identical to observation 0 for all the cells.
        for unit in range(self.n_cells):
            self.sutils_units[unit][1] = self.sutils_units[unit][0]
            
        for ob in range(self.n_observations):
            self.pymuvr_observations.append([])
            for unit in range(self.n_cells):
                self.pymuvr_observations[ob].append(spiketrain_to_list(self.sutils_units[unit][ob]))
        
    def test_compare_with_spykeutils(self):
        sutils_d = stm.van_rossum_multiunit_dist(self.sutils_units,
                                                 weighting=self.cos,
                                                 tau=self.tau)
        pymuvr_d = pymuvr.square_distance_matrix(self.pymuvr_observations,
                                                 self.cos,
                                                 self.tau)
        np.testing.assert_array_almost_equal(sutils_d, pymuvr_d)

if __name__ == "__main__":
    unittest.main(verbosity=2)
