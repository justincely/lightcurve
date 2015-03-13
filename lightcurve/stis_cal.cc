# include <Python.h>
# include <numpy/arrayobject.h>
# include <stdlib.h>

static PyObject *map_image(PyObject *self, PyObject *args) {

  /* local variables */
  int i=0, x=0, y=0;

	PyObject *in_im, *in_x, *in_y;
	PyArrayObject *im, *x_coords, *y_coords;

	if (!PyArg_ParseTuple(args, "OOO", &in_im, &in_x, &in_y)) {
	    PyErr_SetString(PyExc_RuntimeError, "can't read arguments");
	    return NULL;
	}

	x_coords = (PyArrayObject *)PyArray_FROM_OTF(in_x, NPY_INT32, NPY_IN_ARRAY);
	y_coords = (PyArrayObject *)PyArray_FROM_OTF(in_y, NPY_INT32, NPY_IN_ARRAY);
	im = (PyArrayObject *)PyArray_FROM_OTF(in_im, PyArray_TYPE(in_im), NPY_IN_ARRAY);

	if (x_coords == NULL || y_coords == NULL || im == NULL)
	    return NULL;

  int n_rows_in = PyArray_DIM(x_coords, 0);
  npy_intp nrows[1] = {n_rows_in};

  PyArrayObject *out_array;
  out_array = (PyArrayObject *) PyArray_SimpleNew(1, nrows, PyArray_TYPE(in_im));

  for (int i=0; i<n_rows_in; i++){
    x = *(npy_int32 *) PyArray_GETPTR1(x_coords, i);
    y = *(npy_int32 *) PyArray_GETPTR1(y_coords, i);

    if (PyTypeNum_ISFLOAT(PyArray_TYPE(in_im))){
      if ((x >= 2048) || (x <= 0) || (y >= 2048) || (y <= 0)){
        *(npy_float *) PyArray_GETPTR1(out_array, i) = 0.0;
      }
      else{
        *(npy_float *) PyArray_GETPTR1(out_array, i) = *(npy_float *) PyArray_GETPTR2(im, x, y);
      }
    }
    else if (PyTypeNum_ISFLOAT(PyArray_TYPE(in_im))){
      if ((x >= 2048) || (x <= 0) || (y >= 2048) || (y <= 0)){
        *(npy_int *) PyArray_GETPTR1(out_array, i) = 0;
      }
      else{
        *(npy_int *) PyArray_GETPTR1(out_array, i) = *(npy_int *) PyArray_GETPTR2(im, x, y);
      }
    }
  }

  Py_DECREF(x_coords);
  Py_DECREF(y_coords);
  Py_DECREF(im);
  return Py_BuildValue("N", out_array);
}


static PyMethodDef stis_cal_methods[] = {
  {"map_image", map_image, METH_VARARGS, "Map image through coordinates"},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initstis_cal(void) {

	PyObject *mod;		/* the module */
	PyObject *dict;		/* the module's dictionary */

	mod = Py_InitModule("stis_cal", stis_cal_methods);
	import_array();

	/* set the doc string */
	dict = PyModule_GetDict(mod);
	//PyDict_SetItemString(dict, "__doc__",
	//	PyString_FromString(DocString()));
}
