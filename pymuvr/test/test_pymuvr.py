import unittest
import random

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
        mean_isi = 0.1
        max_duration = 3
        self.tau = 1
        self.cos = 0.5
        self.observations = [[simple_train(mean_isi, max_duration) for c in range(n_cells)] for o in range(n_observations)]
        self.observations[0] = self.observations[1][:]
    def test_square_distance_matrix(self):
        d = pymuvr.square_distance_matrix(self.observations, self.cos, self.tau)
        assert d.shape == (len(self.observations), len(self.observations))
        print("")
        print(d)
    def test_distance_matrix(self):
        d = pymuvr.distance_matrix(self.observations[:3],
                                   self.observations[3:],
                                   self.cos,
                                   self.tau)
        assert d.shape == (3, len(self.observations)-3)
        print("")
        print(d)
        

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDistanceMatrix)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
