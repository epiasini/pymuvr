import unittest
import random
from numpy.linalg import norm

import pymuvr

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
        mean_isi = 0.033
        max_duration = 0.5
        self.tau = 0.012
        self.cos = 0.5
        self.observations = [[simple_train(mean_isi, max_duration) for c in range(n_cells)] for o in range(n_observations)]
        self.observations[0] = self.observations[1][:]
    def test_square_distance_matrix(self):
        d = pymuvr.square_distance_matrix(self.observations, self.cos, self.tau)
        assert d.shape == (len(self.observations), len(self.observations))
    def test_distance_matrix(self):
        d = pymuvr.distance_matrix(self.observations[:3],
                                   self.observations[3:],
                                   self.cos,
                                   self.tau)
        assert d.shape == (3, len(self.observations)-3)
    def test_compare_square_and_rectangular(self):
        d_rectangular = pymuvr.distance_matrix(self.observations,
                                               self.observations,
                                               self.cos,
                                               self.tau)
        d_square = pymuvr.square_distance_matrix(self.observations, self.cos, self.tau)
        assert norm(d_rectangular - d_square)/norm((d_rectangular + d_square)/2) < 1e-12

if __name__ == "__main__":
    ## uncomment for verbosity
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestDistanceMatrix)
    #unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()
