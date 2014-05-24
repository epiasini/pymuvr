# pymuvr
## Overview
[![Build Status](https://travis-ci.org/epiasini/pymuvr.svg?branch=master)](https://travis-ci.org/epiasini/pymuvr)

Multi-unit Van Rossum spike train metric. This is a kernel-based
implementation with markage vector and precomputed exponential factor,
as described in Houghton and Kreuz, 2012, 'On the efficient
calculation of Van Rossum distances.'. This package started out as a 
Python wrapping of the original C++ implementation given by the authors
of the paper, and evolved from there with bugfixes and improvements.

## Documentation
Full documentation is hosted at http://pymuvr.readthedocs.org/

## Requirements
pymuvr requires NumPy>=1.7.

## Installation
To install the latest release, run
```shell
pip install pymuvr
```
If you prefer installing from git, use
```shell
git clone https://github.com/epiasini/pymuvr.git
cd pymuvr
python setup.py install
```
Note that you'll get testing and benchbark scripts only if you install
manually (i.e. not via pip).

## Usage
The module exposes two functions:
```python
pymuvr.distance_matrix(parallel_trains_1, parallel_trains_2, cos, tau)
pymuvr.square_distance_matrix(parallel_trains, cos, tau)
```
`distance_matrix` calculates the 'bipartite' (rectangular) dissimilarity
matrix between the trains in `parallel_trains_1` and those in 
`parallel_trains_2`.

`square_distance_matrix` calculates the 'all-to-all' dissimilarity
matrix between each pair of trains in parallel_trains. It's an optimised
form of `distance_matrix(parallel_trains, parallel_trains, cos, tau)`.

They both return their results as a 2D numpy array.

The `parallel_trains` arguments must be thrice-nested lists of spiketimes,
in such a way that `parallel_trains[i][j][k]` represents the time of the
kth spike of the jth cell of the ith train.
`cos` and `tau` are the usual parameters for the multiunit Van Rossum metric.

See `test/test_pymuvr.py` for detailed examples of usage.

## License
This package is licensed under GPLv3. See COPYING for details.

