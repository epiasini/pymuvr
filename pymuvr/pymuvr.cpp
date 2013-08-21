#include "Python.h"
#include <vector>

#include "numpy/arrayobject.h"
#include "Van_Rossum_Multiunit.hpp"

using namespace std;

// prototypes
static PyObject * distance_matrix(PyObject *self, PyObject *args, PyObject *kwds);

// method table
static PyMethodDef pymuvrMethods[] = {
  {"distance_matrix", (PyCFunction)distance_matrix, METH_VARARGS|METH_KEYWORDS, "Distance (dissimilarity) matrix with the Multi-Unit Van Rossum metric."},
  {NULL, NULL, 0, NULL}  /* Sentinel */
};

// module initialisation
PyMODINIT_FUNC
initpymuvr(void) {
  (void) Py_InitModule("pymuvr", pymuvrMethods);
  import_array()
}


// main module function
static PyObject * distance_matrix(PyObject *self, PyObject *args, PyObject *kwds){
  static char *kwlist[] = {"tau", "cos", NULL};
  double tau, cos;
  PyObject *observations;
  Py_ssize_t big_n, big_p;
  PyArrayObject *py_d_matrix;
  double **c_d_matrix;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!dd", kwlist,
				   &PyObject_Type, &observations, &tau, &cos))
    return NULL;
  
  // build the 3D array (observation, cell, spiketime) to be fed to
  // the original metric implementation
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
  npy_intp dims[2] = {big_n, big_n};
  py_d_matrix = (PyArrayObject *)PyArray_SimpleNew(big_n, dims, NPY_DOUBLE);
  int intn = (int)PyArray_DIM(py_d_matrix, 0);
  int intm = (int)PyArray_DIM(py_d_matrix, 1);
  double *a = (double *)PyArray_DATA(py_d_matrix); /* pointer to data as double */
  for (unsigned int i=0; i<intn; i++){
    c_d_matrix[i]=a+i*intm;
  }
  
  // perform the core distance calculations
  d(c_d_matrix, trains, tau, cos);

  // deallocate the memory used by the data structure needed for
  // direct access to the numpy array
  delete [] c_d_matrix;
  return PyArray_Return(py_d_matrix);
}


