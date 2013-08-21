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
        n_observations = 100
        n_cells = 10
        mean_isi = 3
        max_duration = 25
        self.tau = 3
        self.cos = 0.5
        self.observations = [[simple_train(mean_isi, max_duration) for c in range(n_cells)] for o in range(n_observations)]
        self.observations[0] = self.observations[1][:]
    def test_distance_matrix(self):
        d = pymuvr.distance_matrix(self.observations, self.tau, self.cos)
        assert d.shape == (len(self.observations), len(self.observations))
        print("")
        print(d)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDistanceMatrix)
    unittest.TextTestRunner(verbosity=2).run(suite)
   # unittest.main()
