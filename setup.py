from setuptools import setup, Extension
import os
import codecs
import numpy

version = "1.0.1.dev1"

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with codecs.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

ext_module = Extension("pymuvr",
                       sources=[os.path.join("src", "Van_Rossum_Multiunit.cpp"),
                                os.path.join("src", "pymuvr.cpp")],
                       include_dirs=[numpy.get_include(),
                                     "include"])

setup (name = "pymuvr",
       version = version,
       url = 'https://github.com/epiasini/pymuvr',
       description = "Multi-unit Van Rossum spike train metric",
       long_description = long_description,
       install_requires = ['numpy>=1.7'],
       author = "Eugenio Piasini",
       author_email = "e.piasini@ucl.ac.uk",
       license = "GPLv3+",
       ext_modules = [ext_module],
       test_suite = 'tests',
       include_package_data = True)

