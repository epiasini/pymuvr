pymuvr
======

Overview
--------
.. image:: https://travis-ci.org/epiasini/pymuvr.svg?branch=master
    :target: https://travis-ci.org/epiasini/pymuvr
    :alt: build status

A Python package for the fast calculation of Multi-unit Van Rossum
neural spike train metrics, with the kernel-based algorithm described
in Houghton and Kreuz, *On the efficient calculation of Van Rossum
distances* (Network: Computation in Neural Systems, 2012, 23,
48-58). This package started out as a Python wrapping of the original
C++ implementation given by the authors of the paper, and evolved from
there with bugfixes and improvements.

Documentation
-------------
Full documentation is hosted at http://pymuvr.readthedocs.org/.

Requirements
------------
- Python 2.7 or 3.x.
- NumPy>=1.7.
- C++ development tools and Standard Library (package `build-essential` on Debian).
- Python development tools (package `python-dev` on Debian).

Installation
------------
To install the latest release, run::

  pip install pymuvr

If you prefer installing from git, use::

  git clone https://github.com/epiasini/pymuvr.git
  cd pymuvr
  python setup.py install

Usage
-----
The module exposes two functions::

  pymuvr.distance_matrix(observations1, observations2, cos, tau)

::

   pymuvr.square_distance_matrix(observations, cos, tau)

`distance_matrix` calculates the 'bipartite' (rectangular)
dissimilarity matrix between the multi-unit trains in `observations1`
and those in `observations2`.

`square_distance_matrix` calculates the 'all-to-all' dissimilarity
matrix between each pair of trains in parallel_trains. It's an
optimised form of `distance_matrix(observations, observations, cos,
tau)`.

They both return their results as a 2D numpy array.

The `observations` arguments must be thrice-nested lists of
spiketimes, in such a way that `observations[i][j][k]` represents
the time of the kth spike of the jth cell of the ith observation.  `cos` and
`tau` are the usual parameters for the multiunit Van Rossum metric.

See `examples/benchmark_versus_spykeutils.py` for an example of usage
comparing the performance of pymuvr with the pure Python
implementation of the multiunit Van Rossum distance in
`spykeutils <https://github.com/rproepp/spykeutils>`_.

License
-------
This package is licensed under version 3 of the GPL or any later
version. See COPYING for details.

Getting the source
------------------
Source code for pymuvr is hosted at https://github.com/epiasini/pymuvr.
