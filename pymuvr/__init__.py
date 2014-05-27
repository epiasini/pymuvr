__version__ = '(unknown)'
try:
    import pkg_resources
except ImportError:
    # pkg_resources is provided by setuptools. If it's not available,
    # we can't access the package version metadata.
    pass
else:
    try:
        __version__ = pkg_resources.get_distribution('pymuvr').version
    except pkg_resources.DistributionNotFound:
        # This means that pkg_resources is available, but the package
        # is either not installed, or it has been installed manually
        # without setuptools, so we can't access the version metadata.
        pass

# Expose the spike train distance functions defined in the C++
# extension as top-level objects for the package.
from .native.bindings import distance_matrix, square_distance_matrix
