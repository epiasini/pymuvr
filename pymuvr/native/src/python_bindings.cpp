#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <vector>
#include <iostream>

#include "numpy/arrayobject.h"
#include "van_rossum_multiunit.hpp"

using namespace std;

/*==== prototypes ====*/
static PyObject * distance_matrix(PyObject *self, PyObject *args);

static PyObject * square_distance_matrix(PyObject *self, PyObject *args);

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
const char * distance_matrix_docstring = "distance_matrix(observations1, observations2, cos, tau)\n\n\
Return the *bipartite* (rectangular) dissimilarity matrix between the observations in the first and the second list.\n\n\
\
\
:param list observations1,observations2: Two lists of multi-unit spike trains to compare. Each *observations* parameter must be a thrice-nested list of spike times, with `observations[i][j][k]` representing the time of the kth spike of the jth cell of the ith observation.\n\
:param float cos: mixing parameter controlling the interpolation between *labelled-line* mode (cos=0) and *summed-population* mode (cos=1). It corresponds to the cosine of the angle between the vectors used for the euclidean embedding of the multiunit spike trains.\n\
:param float tau: time scale for the exponential kernel, controlling the interpolation between pure *coincidence detection* (tau=0) and *spike count* mode (very large tau). Note that setting tau=0 is always allowed, but there is a range (0, epsilon) of forbidden values that tau is not allowed to assume. The upper bound of this range is proportional to the absolute value of the largest spike time in *observations*, with the proportionality constant being system-dependent. As a rule of thumb tau and the spike times should be within 4 orders of magnitude of each other; for example, if the largest spike time is 10s a value of tau>1ms will be expected. An exception will be raised if tau falls in the forbidden range.\n\
\
:return: A len(observations1) x len(observations2) numpy array containing the distance between each pair of observations that can be formed by taking one observation from *observations1* and one from *observations2*.\n\
:rtype: *numpy.ndarray*\n\
:raises IndexError: if the number of cells in the observations in *observations1* is differnt from that in *observations2*.\n\
:raises OverflowError: if *tau* falls in the forbidden interval.";

const char * square_distance_matrix_docstring = "square_distance_matrix(observations, cos, tau)\n\n\
Return the *all-to-all* (square) dissimilarity matrix for the given list of observations.\n\n\
\
:param list observations: A list of multi-unit spike trains to compare.\n\
:param float cos: mixing parameter controlling the interpolation between *labelled-line* mode (cos=0) and *summed-population* mode (cos=1).\n\
:param float tau: time scale for the exponential kernel, controlling the interpolation between pure *coincidence detection* (tau=0) and *spike count* mode (very large tau).\n\
\
:return: A len(observations) x len(observations) numpy array containing the distance between all possible pairs of observations.\n\
:rtype: *numpy.ndarray*\n\
:raises: **OverflowError** - if *tau* falls in the forbidden interval.\n\
\
Effectively equivalent to *distance_matrix(observations, observations, cos, tau)*, but optimised for speed. See the *distance_matrix* description for details.\n";

/*==== method table ====*/
static PyMethodDef bindings_methods[] = {
  {"distance_matrix", distance_matrix, METH_VARARGS, distance_matrix_docstring},
  {"square_distance_matrix", square_distance_matrix, METH_VARARGS, square_distance_matrix_docstring},
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
    PyObject *module = Py_InitModule("bindings", bindings_methods);
#endif
    import_array()

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("bindings.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

/*==== function implementations ====*/

static PyObject * distance_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  PyObject *observations1, *observations2, *py_d_matrix;
  Py_ssize_t big_n, big_m, big_p;
  npy_intp dims[2];
  PyArray_Descr *descr;
  double **c_d_matrix;

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "OOdd",
			&observations1, &observations2, &cos, &tau))
    return NULL;
  
  /*
    Build the 3D arrays (observation, cell, spiketime) to be fed to
    the original metric implementation. There are big_n observations
    (ie "network activity realisations", or "multi-unit spike trains")
    in the first set and big_m observations in the second set, and
    each of them is composed by big_p single-unit spike trains, but
    each of the single-unit spike train can have a different
    length. This means the 3D arrays don't have a regular shape.
  */
  big_n = PyList_Size(observations1);
  big_m = PyList_Size(observations2);
  big_p = PyList_Size(PyList_GetItem(observations1, (Py_ssize_t)0));
  if (PyList_Size(PyList_GetItem(observations2, (Py_ssize_t)0)) != big_p){
    PyErr_SetString(PyExc_IndexError, "the observations in both lists must have the same number of cells.");
    return NULL;
  }
  vector<vector<vector<double> > > trains1(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > trains2(big_m,vector<vector<double> >(big_p));
  for(Py_ssize_t n=0;n<big_n;++n){
    PyObject *ob = PyList_GetItem(observations1, n);
    for(Py_ssize_t p=0;p<big_p;++p){
      PyObject *cell = PyList_GetItem(ob, p);
      for(Py_ssize_t s=0;s<PyList_Size(cell);++s){
	trains1.at(n).at(p).push_back(PyFloat_AsDouble(PyList_GetItem(cell, s)));
      }
    }
  }
  for(Py_ssize_t n=0;n<big_m;++n){
    PyObject *ob = PyList_GetItem(observations2, n);
    for(Py_ssize_t p=0;p<big_p;++p){
      PyObject *cell = PyList_GetItem(ob, p);
      for(Py_ssize_t s=0;s<PyList_Size(cell);++s){
	trains2.at(n).at(p).push_back(PyFloat_AsDouble(PyList_GetItem(cell, s)));
      }
    }
  }
  
  /*
    Instantiate the dissimilarity (distance) matrix where the results
    will be stored
  */
  dims[0] = big_n;
  dims[1] = big_m;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }
  /* Perform the core distance calculations */
  try{
    d_exp_markage_rect(c_d_matrix, trains1, trains2, tau, cos);
  }
  catch (char const* e){
    PyErr_SetString(PyExc_OverflowError, e);
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

static PyObject * square_distance_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  PyObject *observations, *py_d_matrix;
  Py_ssize_t big_n, big_p;
  npy_intp dims[2];
  PyArray_Descr *descr;
  double **c_d_matrix;

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "Odd",
			&observations, &cos, &tau))
    return NULL;
  
  /*
    Build the 3D array (observation, cell, spiketime) to be fed to the
    original metric implementation. There are big_n observations (ie
    "network activity realisations", or "multi-unit spike trains"),
    and each of them is composed by big_p single-unit spike trains,
    but each of the single-unit spike train can have a different
    length. This means the 3D array doesn't have a regular shape.
  */
  big_n = PyList_Size(observations);
  big_p = PyList_Size(PyList_GetItem(observations, (Py_ssize_t)0));
  vector<vector<vector<double> > > trains(big_n,vector<vector<double> >(big_p));
  for(Py_ssize_t n=0;n<big_n;++n){
    PyObject *ob = PyList_GetItem(observations, n);
    for(Py_ssize_t p=0;p<big_p;++p){
      PyObject *cell = PyList_GetItem(ob, p);
      for(Py_ssize_t s=0;s<PyList_Size(cell);++s){
	trains.at(n).at(p).push_back(PyFloat_AsDouble(PyList_GetItem(cell, s)));
      }
    }
  }
  
  /*
    Instantiate the dissimilarity (distance) matrix where the results
    will be stored
  */
  dims[0] = big_n;
  dims[1] = big_n;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }

  /* Perform the core distance calculations */
  try{
  d_exp_markage(c_d_matrix, trains, tau, cos);
  }
  catch (char const* e){
    PyErr_SetString(PyExc_OverflowError, e);
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


