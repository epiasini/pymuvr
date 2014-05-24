Usage
=====

The module exposes two functions::

  pymuvr.distance_matrix(parallel_trains_1, parallel_trains_2, cos, tau)

::

  pymuvr.square_distance_matrix(parallel_trains, cos, tau)

``distance_matrix`` calculates the *bipartite* (rectangular) dissimilarity
matrix between the trains in ``parallel_trains_1`` and those in 
``parallel_trains_2``.

``square_distance_matrix`` calculates the *all-to-all* dissimilarity
matrix between each pair of trains in parallel_trains. It's an optimised
form of ``distance_matrix(parallel_trains, parallel_trains, cos, tau)``.

They both return their results as a 2D numpy array.

The ``parallel_trains`` arguments must be thrice-nested lists of spiketimes,
in such a way that ``parallel_trains[i][j][k]`` represents the time of the
kth spike of the jth cell of the ith train.
``cos`` and ``tau`` are the usual parameters for the multiunit Van Rossum metric.

See ``test/test_pymuvr.py`` for detailed examples of usage.
