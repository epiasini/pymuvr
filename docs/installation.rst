Installation
============

Requirements
------------

- CPython 2.7 or 3.x.
- NumPy>=1.7.
- C++ development tools and Standard Library (package `build-essential` on Debian).
- Python development tools (package `python-dev` on Debian).

Installing via *pip*
--------------------
To install the latest release, run::

  pip install pymuvr

Testing
-------

From the root directory of the source distribution, run::

  python setup.py test

(requires *setuptools*). **Alternatively**, if pymuvr is already
installed on your system, look for the copy of the ``test_pymuvr.py``
script installed alongside the rest of the pymuvr files and execute
it. For example::

  python /usr/lib/pythonX.Y/site-packages/pymuvr/test/test_pymuvr.py
