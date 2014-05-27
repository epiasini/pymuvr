#!/bin/bash

# This script can be used for testing with Travis CI. If the
# with_spikeutils variable is true, we try to install spykeutils from
# pip. This could still fail, but the tests should pass anyway as
# pymuvr does not depend on spykeutils. Otherwise, we just run the
# test suite without even attempting to install spykeutils. This
# allows for automatically testing with and without spykeutils on
# platforms where it is supported, while allowing for the installation
# to fail where it's not without erroring the test setup.

if [ "$with_spykeutils" = true ]; then
    echo "trying to install spykeutils"
    # spykeutils depends on scipy, which depends on libatlas.
    sudo apt-get update -qq
    sudo apt-get install -qq liblapack-dev libatlas-dev gfortran
    
    # Install spykeutils. This could either succeed or fail depending
    # on the Python version, as at the time of writing spykeutils is
    # not compatible with Python 3. Anyway, the test script will just
    # skip the relevant tests if it finds it can't import spykeutils.
    pip install -q spykeutils
else
    echo "not trying to install spykeutils"
fi

python pymuvr/test/test_pymuvr.py
