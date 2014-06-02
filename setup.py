try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension
import os
import codecs
import numpy

version = '1.2.0'

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with codecs.open(os.path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

ext_module = Extension("pymuvr.native.bindings",
                       sources=["pymuvr/native/src/van_rossum_multiunit.cpp",
                                "pymuvr/native/src/python_bindings.cpp"],
                       include_dirs=[numpy.get_include(),
                                     "pymuvr/native/include"])

setup (name="pymuvr",
       version=version,
       url="https://github.com/epiasini/pymuvr",
       description="Multi-unit Van Rossum spike train metric",
       long_description=long_description,
       install_requires=["numpy>=1.7"],
       author="Eugenio Piasini",
       author_email="e.piasini@ucl.ac.uk",
       license="GPLv3+",
       classifiers=[
           "Development Status :: 5 - Production/Stable",
           "Intended Audience :: Science/Research",
           "Topic :: Scientific/Engineering",
           "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
           "Programming Language :: Python :: Implementation :: CPython",
           "Programming Language :: Python :: 2.7",
           "Programming Language :: Python :: 3.2",
           "Programming Language :: Python :: 3.3"
       ],
       packages=["pymuvr",
                 "pymuvr.native",
                 "pymuvr.test"],
       ext_modules=[ext_module],
       # examples must be included as package_data
       package_data={"pymuvr": ["examples/*.py"]},
       test_suite="pymuvr.test",
       zip_safe=False)

