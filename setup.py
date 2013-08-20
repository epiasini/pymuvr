from distutils.core import setup, Extension

module1 = Extension('pymuvr.pymuvr',
                    sources = ['pymuvr/Van_Rossum_Multiunit.cpp', 'pymuvr/pymuvr.cpp'],
                    include_dirs = ['include'])

setup (name = 'pymuvr',
       version = '0.1',
       description = 'Multi-unit Van Rossum spike train metric',
       ext_modules = [module1],
       long_description = "Multi-unit Van Rossum spike train metric, kernel-based version. This is a Python wrapping of the C++ implementation given by the authors in Houghton and Kreuz, 2012, 'On the efficient calculation of Van Rossum distances.'")
