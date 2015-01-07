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

Full documentation (with examples and full reference) is hosted at
http://pymuvr.readthedocs.org/.

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

License
-------
This package is licensed under version 3 of the GPL or any later
version. See COPYING for details.

Getting the source
------------------
Source code for pymuvr is hosted at https://github.com/epiasini/pymuvr.
