from distutils.core import setup, Extension
import numpy

module1 = Extension("pymuvr",
                    sources = ["src/Van_Rossum_Multiunit.cpp",
                               "src/pymuvr.cpp"],
                    include_dirs = [numpy.get_include(),
                                    "include"])

setup (name = "pymuvr",
       version = "1.0.1",
       url = "https://github.com/epiasini/pymuvr",
       description = "Multi-unit Van Rossum spike train metric",
       long_description = "Multi-unit Van Rossum spike train metric; kernel-based version with markage vector and precomputed exponential factor, as described in Houghton and Kreuz, 2012, 'On the efficient calculation of Van Rossum distances.'. This is a Python wrapping of the original C++ implementation given by the authors of the paper.",
       install_requires = ['numpy>=1.7'],
       author = "Eugenio Piasini",
       author_email = "e.piasini@ucl.ac.uk",
       license = "GPLv3 or later",
       ext_modules = [module1])

