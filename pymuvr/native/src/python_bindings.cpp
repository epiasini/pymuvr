#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <vector>
#include <stdexcept>

#include "numpy/arrayobject.h"
#include "convolved_spike_train.hpp"
#include "van_rossum_multiunit.hpp"

using std::vector; //from vector
using std::overflow_error; using std::invalid_argument; // from stdexcept

/*==== prototypes ====*/
static PyObject * dissimilarity_matrix(PyObject *self, PyObject *args);

static PyObject * square_dissimilarity_matrix(PyObject *self, PyObject *args);

/*!
 * Convert spike trains from their Python-based representation to a
 * data structure based on ConvolvedSpikeTrain.
 *
 * At the Python level, sets of multiunit spike trains are represented
 * as thrice-nested lists, where observations[i][j][k] is the kth
 * spiketime of the jth cell in the ith observation. Internally, we
 * use a vector-of-vector of instances of ConvolvedSpikeTrain. This
 * function converts between the two data structures by extracting
 * spike time data from the Python list-of-lists-of-lists and using it
 * to fill the vector-of-vectors out-parameter with instances of
 * ConvolvedSpikeTrain.
 *
 * \param[in] observations A Python list-of-lists-of-lists object
 * representing a set of multiunit spike trains, such that \a
 * observations[i][j][k] is the kth spiketime of the jth cell in the
 * ith observation.
 *
 * \param[in] tau Time scale for the exponential kernel.
 *
 * \param[out] trains 2D array (observation, cell) that will be filled
 * with convolved spike trains.
 */
static void nested_lists_to_spiketrain_vectors(PyObject *observations, double tau, vector <vector <ConvolvedSpikeTrain> > &trains);


/*==== module initialisation ====*/
/*
  See https://wiki.python.org/moin/PortingExtensionModulesToPy3k for
  the template used to support Python 2 and Python 3 at the same time.
*/
struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

/*==== docstrings ====*/
const char * dissimilarity_matrix_docstring = "dissimilarity_matrix(observations1, observations2, cos, tau, mode)\n\n\
Return the *bipartite* (rectangular) dissimilarity matrix between the observations in the first and the second list.\n\n\
\
\
:param list observations1,observations2: Two lists of multi-unit spike trains to compare. Each *observations* parameter must be a thrice-nested list of spike times, with `observations[i][j][k]` representing the time of the kth spike of the jth cell of the ith observation.\n\
:param float cos: mixing parameter controlling the interpolation between *labelled-line* mode (cos=0) and *summed-population* mode (cos=1). It corresponds to the cosine of the angle between the vectors used for the euclidean embedding of the multiunit spike trains.\n\
:param float tau: time scale for the exponential kernel, controlling the interpolation between pure *coincidence detection* (tau=0) and *spike count* mode (very large tau). Note that setting tau=0 is always allowed, but there is a range (0, epsilon) of forbidden values that tau is not allowed to assume. The upper bound of this range is proportional to the absolute value of the largest spike time in *observations*, with the proportionality constant being system-dependent. As a rule of thumb tau and the spike times should be within 4 orders of magnitude of each other; for example, if the largest spike time is 10s a value of tau>1ms will be expected. An exception will be raised if tau falls in the forbidden range.\n\
\
:param string mode: type of dissimilarity measure to be computed. Must be either 'distance' or 'inner product'.\n\
:return: A len(observations1) x len(observations2) numpy array containing the dissimilarity (distance or inner product) between each pair of observations that can be formed by taking one observation from *observations1* and one from *observations2*.\n\
:rtype: *numpy.ndarray*\n\
:raises IndexError: if the observations in *observations1* and *observations2* don't have all the same number of cells.\n\
:raises OverflowError: if *tau* falls in the forbidden interval.";

const char * square_dissimilarity_matrix_docstring = "square_dissimilarity_matrix(observations, cos, tau, mode)\n\n\
Return the *all-to-all* (square) dissimilarity matrix for the given list of observations.\n\n\
\
:param list observations: A list of multi-unit spike trains to compare.\n\
:param float cos: mixing parameter controlling the interpolation between *labelled-line* mode (cos=0) and *summed-population* mode (cos=1).\n\
:param float tau: time scale for the exponential kernel, controlling the interpolation between pure *coincidence detection* (tau=0) and *spike count* mode (very large tau).\n\
\
:param string mode: type of dissimilarity measure to be computed. Must be either 'distance' or 'inner product'.\n\
:return: A len(observations) x len(observations) numpy array containing the dissimilarity (distance or inner product) between all possible pairs of observations.\n\
:rtype: *numpy.ndarray*\n\
:raises IndexError: if the observations in *observations* don't have all the same number of cells.\n\
:raises OverflowError: if *tau* falls in the forbidden interval.\n\
\
Effectively equivalent to ``dissimilarity_matrix(observations, observations, cos, tau)``, but optimised for speed. See :func:`pymuvr.dissimilarity_matrix` for details.";

/*==== method table ====*/
static PyMethodDef bindings_methods[] = {
  {"dissimilarity_matrix", dissimilarity_matrix, METH_VARARGS, dissimilarity_matrix_docstring},
  {"square_dissimilarity_matrix", square_dissimilarity_matrix, METH_VARARGS, square_dissimilarity_matrix_docstring},
  {NULL, NULL, 0, NULL}  // Sentinel - marks the end of the structure
};

/*==== back to module initialisation ====*/
#if PY_MAJOR_VERSION >= 3

static int bindings_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int bindings_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "bindings",
        NULL,
        sizeof(struct module_state),
        bindings_methods,
        NULL,
        bindings_traverse,
        bindings_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_bindings(void)

#else
#define INITERROR return

PyMODINIT_FUNC
initbindings(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule((char *) "bindings", bindings_methods);
#endif
    import_array()

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException((char *) "bindings.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

/*==== function implementations ====*/

static PyObject * dissimilarity_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  char *mode;
  PyObject *observations1, *observations2, *py_d_matrix;
  Py_ssize_t big_n, big_m, big_p;
  npy_intp dims[2];
  PyArray_Descr *descr;
  double **c_d_matrix;

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "OOdds",
			&observations1, &observations2, &cos, &tau, &mode))
    return NULL;
  
  /* Find out the number of observations per set and the number of
     cells per observation. */
  big_n = PyList_Size(observations1);
  big_m = PyList_Size(observations2);
  big_p = PyList_Size(PyList_GetItem(observations1, (Py_ssize_t)0));

  /* Check that number of cells in observations1 matches that in observations2 */
  if (PyList_Size(PyList_GetItem(observations2, (Py_ssize_t)0)) != big_p){
    PyErr_SetString(PyExc_IndexError, "trying to compare observations with a different number of cells.");
    return NULL;
  }

  /*
    Build the 2D arrays (observation, cell) of convolved spike trains
    to be fed to the dissimilarity function.
  */
  vector<vector <ConvolvedSpikeTrain> > trains1(big_n, vector<ConvolvedSpikeTrain>(big_p));
  vector<vector <ConvolvedSpikeTrain> > trains2(big_m, vector<ConvolvedSpikeTrain>(big_p));
  nested_lists_to_spiketrain_vectors(observations1, tau, trains1);
  nested_lists_to_spiketrain_vectors(observations2, tau, trains2);
  /* Check if something went wrong while generating the convolved spike trains */
  if (PyErr_Occurred() != NULL) {
    return NULL;
  }
  
  /*
    Instantiate the dissimilarity matrix where the results will be
    stored
  */
  dims[0] = big_n;
  dims[1] = big_m;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }
  /* Perform the core dissimilarity calculations */
  try{
    if (strcmp(mode, "inner product")==0) {
      inner_product(trains1, trains2, cos, c_d_matrix);
    } else if (strcmp(mode, "distance")==0) {
      distance(trains1, trains2, cos, c_d_matrix);
    } else {
      PyErr_SetString(PyExc_ValueError, "please specify the type of dissimilarity function ('distance' or 'inner product') to calculate.");
      goto fail;
    }
  }catch (invalid_argument const& e){
    PyErr_SetString(PyExc_IndexError, e.what());
    goto fail;
  }

  /*
    Free the memory used by the data structure needed for direct
    access to the numpy array
  */
  PyArray_Free(py_d_matrix, c_d_matrix);

  return PyArray_Return((PyArrayObject *)py_d_matrix);

 fail:
  /*
    If something goes wrong, remember to free the memory used for
    c_d_matrix and to decrease the reference count for py_d_matrix
    before returning.
  */
  PyArray_Free(py_d_matrix, c_d_matrix);
  Py_DECREF(py_d_matrix);
  return NULL;
}

static PyObject * square_dissimilarity_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  char *mode;
  PyObject *observations, *py_d_matrix;
  Py_ssize_t big_n, big_p;
  npy_intp dims[2];
  PyArray_Descr *descr;
  double **c_d_matrix;

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "Odds",
			&observations, &cos, &tau, &mode))
    return NULL;
  
  /* Find out the number of observations and the number of cells per
     observation. */
  big_n = PyList_Size(observations);
  big_p = PyList_Size(PyList_GetItem(observations, (Py_ssize_t)0));

  /*
    Build the 2D array (observation, cell) of convolved spike trains
    to be fed to the dissimilarity function.
  */
  vector<vector <ConvolvedSpikeTrain> > trains(big_n, vector<ConvolvedSpikeTrain>(big_p));
  nested_lists_to_spiketrain_vectors(observations, tau, trains);
  /* Check if something went wrong while generating the convolved spike trains */
  if (PyErr_Occurred() != NULL) {
    return NULL;
  }
  
  /*
    Instantiate the dissimilarity matrix where the results will be
    stored
  */
  dims[0] = big_n;
  dims[1] = big_n;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }

  /* Perform the core dissimilarity calculations */
  try{
    if (strcmp(mode, "inner product")==0) {
      inner_product(trains, cos, c_d_matrix);
    } else if (strcmp(mode, "distance")==0) {
      distance(trains, cos, c_d_matrix);
    } else {
      PyErr_SetString(PyExc_ValueError, "please specify the type of dissimilarity function ('distance' or 'inner product') to calculate.");
      goto fail;
    }
  }catch (invalid_argument const& e){
    PyErr_SetString(PyExc_ValueError, e.what());
    goto fail;
  }
  /*
    Free the memory used by the data structure needed for direct
    access to the numpy array
  */
  PyArray_Free(py_d_matrix, c_d_matrix);
  return PyArray_Return((PyArrayObject *)py_d_matrix);

 fail:
  /*
    If something goes wrong, remember to free the memory used for
    c_d_matrix and to decrease the reference count for py_d_matrix
    before returning.
  */
  PyArray_Free(py_d_matrix, c_d_matrix);
  Py_DECREF(py_d_matrix);
  return NULL;
}


static void nested_lists_to_spiketrain_vectors(PyObject *observations, double tau, vector <vector <ConvolvedSpikeTrain> > &trains)
{
  Py_ssize_t big_n = PyList_Size(observations);
  Py_ssize_t big_p = PyList_Size(PyList_GetItem(observations, (Py_ssize_t)0));
  /* Check that input observations all have the same number of
     cells. There are big_n observations (ie "network activity
     realisations", or "multi-unit spike trains"), and each of them
     must be composed composed by big_p single-unit spike trains, but
     each of the single-unit spike train can have a different
     length. */
  for(Py_ssize_t n=0;n<big_n;++n){
    if (PyList_Size(PyList_GetItem(observations, n)) != big_p){
      PyErr_SetString(PyExc_IndexError, "trying to compare observations with a different number of cells.");
    }
  }

  /* Copy data from Python structure to our internal representation */
  for(Py_ssize_t n=0;n<big_n;++n){
    PyObject *ob = PyList_GetItem(observations, n);
    for(Py_ssize_t p=0;p<big_p;++p){
      PyObject *cell = PyList_GetItem(ob, p);
      vector <double> spikes; 
      for(Py_ssize_t s=0;s<PyList_Size(cell);++s){
	spikes.push_back(PyFloat_AsDouble(PyList_GetItem(cell, s)));
      }
      try {
	trains.at(n).at(p) = ConvolvedSpikeTrain(spikes, tau);
      } catch (overflow_error const& e) {
	PyErr_SetString(PyExc_OverflowError, e.what());
      }
    }
  }
}
