#include "Python.h"
#include <vector>

#include "numpy/arrayobject.h"
#include "Van_Rossum_Multiunit.hpp"

using namespace std;

// prototypes
static PyObject * distance_matrix(PyObject *self, PyObject *args);

// method table
static PyMethodDef pymuvrMethods[] = {
  {"distance_matrix", distance_matrix, METH_VARARGS, "Distance (dissimilarity) matrix with the Multi-Unit Van Rossum metric."},
  {NULL, NULL, 0, NULL}  /* Sentinel */
};

// module initialisation
PyMODINIT_FUNC
initpymuvr(void) {
  (void) Py_InitModule("pymuvr", pymuvrMethods);
  import_array()
}


// main module function
static PyObject * distance_matrix(PyObject *self, PyObject *args){
  double cos, tau;
  PyObject *observations;
  Py_ssize_t big_n, big_p;
  PyArrayObject *py_d_matrix;
  double *py_data;
  npy_intp dims[2];
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
  py_d_matrix = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  py_data = (double *)PyArray_DATA(py_d_matrix); // pointer to data as double
  c_d_matrix = new double* [big_n]; // array of pointers to the rows
				    // of the numpy matrix
  for (npy_intp i=0; i<big_n; i++){
    c_d_matrix[i] = py_data + i*big_n;
  }
  
  // perform the core distance calculations
  d_exp_markage(c_d_matrix, trains, tau, cos);

  // deallocate the memory used by the data structure needed for
  // direct access to the numpy array
  delete [] c_d_matrix;

  
  return PyArray_Return(py_d_matrix);
}


