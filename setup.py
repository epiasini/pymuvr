from distutils.core import setup, Extension

module1 = Extension("pymuvr",
                    sources = ["pymuvr/Van_Rossum_Multiunit.cpp", "pymuvr/pymuvr.cpp"])

setup (name = "pymuvr",
       version = "1.0",
       description = "Multi-unit Van Rossum spike train metric",
       long_description = "Multi-unit Van Rossum spike train metric, kernel-based version with markage vector and precomputed exponential factor, as described in Houghton and Kreuz, 2012, 'On the efficient calculation of Van Rossum distances.'. This is a Python wrapping of the original C++ implementation given by the authors of the paper.",
       author = "Eugenio Piasini",
       author_email = "e.piasini@ucl.ac.uk",
       license = "GPLv3 or later",
       ext_modules = [module1])

