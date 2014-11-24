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
from .native.bindings import dissimilarity_matrix, square_dissimilarity_matrix

def distance_matrix(trains1, trains2, cos, tau):
    """
    Return the *bipartite* (rectangular) distance matrix between the observations in the first and the second list.

    Convenience function; equivalent to ``dissimilarity_matrix(trains1, trains2, cos, tau, "distance")``. Refer to :func:`pymuvr.dissimilarity_matrix` for full documentation.
    """
    return dissimilarity_matrix(trains1, trains2, cos, tau, "distance")

def square_distance_matrix(trains, cos, tau):
    """
    Return the *all-to-all* (square) distance matrix for the given list of observations.

    Convenience function; equivalent to ``square_dissimilarity_matrix(trains, cos, tau, "distance")``. Refer to :func:`pymuvr.square_dissimilarity_matrix` for full documentation.
    """
    return square_dissimilarity_matrix(trains, cos, tau, "distance")
