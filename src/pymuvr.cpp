#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include <vector>
#include <iostream>

#include "numpy/arrayobject.h"
#include "Van_Rossum_Multiunit.hpp"

using namespace std;

/*==== prototypes ====*/
static PyObject * distance_matrix(PyObject *self, PyObject *args);

static PyObject * square_distance_matrix(PyObject *self, PyObject *args);

/*==== method table ====*/
static PyMethodDef pymuvrMethods[] = {
  {"distance_matrix", distance_matrix, METH_VARARGS, "Return the 'bipartite' rectangular dissimilarity matrix between the observations in the first and the second list."},
  {"square_distance_matrix", square_distance_matrix, METH_VARARGS, "Return the all-to-all dissimilarity matrix for the given list of observations."},
  {NULL, NULL, 0, NULL}  // Sentinel - marks the end of the structure
};

/*==== module initialisation ====*/
PyMODINIT_FUNC
initpymuvr(void) {
  (void) Py_InitModule("pymuvr", pymuvrMethods);
  import_array();
}


/*==== function implementations ====*/

static PyObject * distance_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  PyObject *observations1, *observations2, *py_d_matrix;
  Py_ssize_t big_n, big_m, big_p;
  npy_intp dims[2];
  PyArray_Descr *descr;
  double **c_d_matrix;

  // parse arguments
  if (!PyArg_ParseTuple(args, "OOdd",
			&observations1, &observations2, &cos, &tau))
    return NULL;
  
  // build the 3D arrays (observation, cell, spiketime) to be fed to
  // the original metric implementation. There are big_n observations
  // (ie "network activity realisations", or "multi-unit spike
  // trains") in the first set and big_m observations in the second
  // set, and each of them is composed by big_p single-unit spike
  // trains, but each of the single-unit spike train can have a
  // different length. This means the 3D arrays don't have a regular
  // shape.
  big_n = PyList_Size(observations1);
  big_m = PyList_Size(observations2);
  big_p = PyList_Size(PyList_GetItem(observations1, (Py_ssize_t)0));
  if (PyList_Size(PyList_GetItem(observations2, (Py_ssize_t)0)) != big_p){
    PyErr_SetString(PyExc_RuntimeError, "pymuvr: the observations in both lists must have the same number of cells.");
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
  
  // instantiate the dissimilarity (distance) matrix where the results
  // will be stored 
  dims[0] = big_n;
  dims[1] = big_m;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "pymuvr: something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }
  // perform the core distance calculations
  try{
    d_exp_markage_rect(c_d_matrix, trains1, trains2, tau, cos);
  }
  catch (char const* e){
    PyErr_SetString(PyExc_RuntimeError, e);
    goto fail;
  }

  // free the memory used by the data structure needed for direct
  // access to the numpy array
  PyArray_Free(py_d_matrix, c_d_matrix);

  return PyArray_Return((PyArrayObject *)py_d_matrix);

 fail:
  // if something goes wrong, remember to free the memory used for
  // c_d_matrix and to decrease the reference count for py_d_matrix
  // before returning.
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

  // parse arguments
  if (!PyArg_ParseTuple(args, "Odd",
			&observations, &cos, &tau))
    return NULL;
  
  // build the 3D array (observation, cell, spiketime) to be fed to
  // the original metric implementation. There are big_n observations
  // (ie "network activity realisations", or "multi-unit spike
  // trains"), and each of them is composed by big_p single-unit spike
  // trains, but each of the single-unit spike train can have a
  // different length. This means the 3D array doesn't have a regular
  // shape.
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
  
  // instantiate the dissimilarity (distance) matrix where the results
  // will be stored 
  dims[0] = big_n;
  dims[1] = big_n;
  py_d_matrix = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  descr = PyArray_DescrFromType(NPY_DOUBLE);
  if (PyArray_AsCArray (&py_d_matrix, &c_d_matrix, dims, 2, descr) == -1){
    PyErr_SetString(PyExc_RuntimeError, "pymuvr: something went wrong while allocating the dissimilarity matrix.");
    goto fail;
  }

  // perform the core distance calculations
  try{
  d_exp_markage(c_d_matrix, trains, tau, cos);
  }
  catch (char const* e){
    PyErr_SetString(PyExc_RuntimeError, e);
    goto fail;
  }
  // free the memory used by the data structure needed for direct
  // access to the numpy array
  PyArray_Free(py_d_matrix, c_d_matrix);
  return PyArray_Return((PyArrayObject *)py_d_matrix);

 fail:
  // if something goes wrong, remember to free the memory used for
  // c_d_matrix and to decrease the reference count for py_d_matrix
  // before returning.
  PyArray_Free(py_d_matrix, c_d_matrix);
  Py_DECREF(py_d_matrix);
  return NULL;
}


