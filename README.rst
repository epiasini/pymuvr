pymuvr
======

Overview
--------
.. image:: https://travis-ci.org/epiasini/pymuvr.svg?branch=master
    :target: https://travis-ci.org/epiasini/pymuvr
    :alt: build status

Multi-unit Van Rossum spike train metric. This is a kernel-based
implementation with markage vector and precomputed exponential factor,
as described in Houghton and Kreuz, 2012, 'On the efficient
calculation of Van Rossum distances.'. This package started out as a
Python wrapping of the original C++ implementation given by the
authors of the paper, and evolved from there with bugfixes and
improvements.

Documentation
-------------
Full documentation is hosted at http://pymuvr.readthedocs.org/

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

Note that you'll get testing and benchbark scripts only if you install
manually (i.e. not via pip).

Usage
-----
The module exposes two functions::

  pymuvr.distance_matrix(parallel_trains_1, parallel_trains_2, cos, tau)

::

   pymuvr.square_distance_matrix(parallel_trains, cos, tau)

`distance_matrix` calculates the 'bipartite' (rectangular)
dissimilarity matrix between the trains in `parallel_trains_1` and
those in `parallel_trains_2`.

`square_distance_matrix` calculates the 'all-to-all' dissimilarity
matrix between each pair of trains in parallel_trains. It's an
optimised form of `distance_matrix(parallel_trains, parallel_trains,
cos, tau)`.

They both return their results as a 2D numpy array.

The `parallel_trains` arguments must be thrice-nested lists of
spiketimes, in such a way that `parallel_trains[i][j][k]` represents
the time of the kth spike of the jth cell of the ith train.  `cos` and
`tau` are the usual parameters for the multiunit Van Rossum metric.

See `examples/benchmark_versus_spykeutils.py` for an example of usage
comparing the performance of pymuvr with the pure Python
implementation of the multiunit Van Rossum distance in
`spykeutils <https://github.com/rproepp/spykeutils>`_.

License
-------
This package is licensed under version 3 of the GPL or any later
version. See COPYING for details.

