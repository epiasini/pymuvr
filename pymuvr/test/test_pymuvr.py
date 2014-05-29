import unittest
import random
import os
import numpy as np

import pymuvr

# check if spykeutils is available. As a special case needed for
# Travis, pretend it's not available in any case if the
# without_spykeutils environment variable is set and is not an empty
# string.
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

class TestTrivialTrains(unittest.TestCase):
    def setUp(self):
        self.tau = 0.012
        self.cos = 0.5

    def test_empty_spike_trains(self):
        observations = [[[]], [[]]]
        d = pymuvr.distance_matrix(observations,
                                   observations,
                                   self.cos, self.tau)
        np.testing.assert_array_equal(d, np.zeros_like(d))

    def test_identical_trains(self):
        observations = [[[1.,2.],[1.5]], [[1.,2.],[1.5]]]
        d = pymuvr.distance_matrix(observations,
                                   observations,
                                   self.cos, self.tau)
        np.testing.assert_array_equal(d, np.zeros_like(d))

    def test_missing_spike(self):
        observations = [[[1.,2.]], [[1.]]]
        d = pymuvr.distance_matrix(observations,
                                   observations,
                                   self.cos, self.tau)
        np.testing.assert_array_equal(d, np.array([[0,1],[1,0]]))

    def test_small_tau_limit(self):
        observations = [[[0.1, 0.2, 0.3]],
                        [[0.1, 0.2, 0.3]],
                        [[0.1, 0.2, 0.4]],
                        [[0.1, 0.2]]]
        d_small = pymuvr.square_distance_matrix(observations, self.cos, 1e-4)
        d_zero = pymuvr.square_distance_matrix(observations, self.cos, 0)
        np.testing.assert_allclose(d_small, d_zero)

    def test_zero_tau(self):
        observations = [[[0.1, 0.2, 0.3]],
                        [[0.1, 0.2, 0.3]],
                        [[0.1, 0.2, 0.4]],
                        [[0.1, 0.2]]]
        target_d = np.array([[0, 0, np.sqrt(2), 1],
                             [0, 0, np.sqrt(2), 1],
                             [np.sqrt(2), np.sqrt(2), 0, 1],
                             [1, 1, 1, 0]])
        d = pymuvr.distance_matrix(observations, observations, 0, 0)
        np.testing.assert_array_equal(d, target_d)

    def test_large_tau_limit(self):
        observations = [[[1,2]], [[1,2]], [[1]]]
        d = pymuvr.square_distance_matrix(observations, self.cos, 1e20)
        np.testing.assert_allclose(d, np.array([[0,0,1],[0,0,1],[1,1,0]]))
        
class TestRandomTrains(unittest.TestCase):
    def setUp(self):
        self.n_observations = 10
        self.n_cells = 20
        self.mean_isi = 0.03
        self.max_duration = 2
        self.cos = np.linspace(0, 1, 3)
        self.tau = np.linspace(0, 0.018, 3)
        self.observations = [[simple_train(self.mean_isi, self.max_duration) for c in range(self.n_cells)] for o in range(self.n_observations)]
        # observation 1 is identical to observation 0 for all the cells.
        self.observations[0] = self.observations[1][:]

    def test_square_distance_matrix(self):
        for cos in self.cos:
            for tau in self.tau:
                d = pymuvr.square_distance_matrix(self.observations, cos, tau)
                self.assertEqual(d.shape, (len(self.observations), len(self.observations)))

    def test_distance_matrix(self):
        for cos in self.cos:
            for tau in self.tau:
                d = pymuvr.distance_matrix(self.observations[:3],
                                           self.observations[3:],
                                           cos,
                                           tau)
                self.assertEqual(d.shape, (3, len(self.observations)-3))

    def test_compare_square_and_rectangular(self):
        for cos in self.cos:
            for tau in self.tau:
                d_rectangular = pymuvr.distance_matrix(self.observations,
                                                       self.observations,
                                                       cos,
                                                       tau)
                d_square = pymuvr.square_distance_matrix(self.observations,
                                                         cos,
                                                         tau)
                np.testing.assert_allclose(d_rectangular, d_square, atol=5e-5)

    def test_empty_spike_train(self):
        observations = [o[:] for o in self.observations]
        observations[0][0] = []
        for cos in self.cos:
            for tau in self.tau:
                d_rectangular = pymuvr.distance_matrix(observations[:3],
                                                       observations[3:],
                                                       cos,
                                                       tau)


@unittest.skipIf(not SPYKEUTILS_IS_AVAILABLE,
                 "can't import spykeutils")
class TestCompareWithSpykeutils(unittest.TestCase):
    def setUp(self):
        self.n_observations = 10
        self.n_cells = 5
        self.rate = 30
        self.tstop = 1
        self.cos = np.linspace(0, 1, 5)
        self.tau = np.linspace(0.0001, 0.018, 3)
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
                self.pymuvr_observations[ob].append(self.sutils_units[unit][ob].tolist())
        
    def test_compare_rectangular_with_spykeutils(self):
        for cos in self.cos:
            for tau in self.tau:
                sutils_d = stm.van_rossum_multiunit_dist(self.sutils_units,
                                                         weighting=cos,
                                                         tau=tau)
                pymuvr_d = pymuvr.distance_matrix(self.pymuvr_observations,
                                                  self.pymuvr_observations,
                                                  cos,
                                                  tau)
                np.testing.assert_allclose(sutils_d, pymuvr_d, atol=5e-5)

    def test_compare_square_with_spykeutils(self):
        for cos in self.cos:
            for tau in self.tau:
                sutils_d = stm.van_rossum_multiunit_dist(self.sutils_units,
                                                         weighting=cos,
                                                         tau=tau)
                pymuvr_d_square = pymuvr.square_distance_matrix(self.pymuvr_observations,
                                                                cos,
                                                                tau)                
                np.testing.assert_allclose(sutils_d, pymuvr_d_square, atol=5e-5)

if __name__ == "__main__":
    unittest.main(verbosity=2)
