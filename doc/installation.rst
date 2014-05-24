Installation
============

Requirements
------------

- Python 2.7 or 3.x.
- NumPy>=1.7.
- C++ development tools and Standard Library (package `build-essential` on Debian).
- Python development tools (package `python-dev` on Debian).

Using *pip*
-----------
To install the latest release, run::

  pip install pymuvr

Using *git* and *setuptools*
----------------------------
If you prefer installing from git, use::

  git clone https://github.com/epiasini/pymuvr.git
  cd pymuvr
  python setup.py install

Note that you'll get testing and benchbark scripts only if you install
manually (i.e. not via pip).
