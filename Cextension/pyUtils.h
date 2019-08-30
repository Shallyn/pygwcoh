/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_PYUTILS__
#define __INCLUDE_PYUTILS__

#include <Python.h>
#include <numpy/arrayobject.h>
#include "cohUtils.h"

#define MACRO_PyArray2REAL8Vector(py_obj, vec) _PyArray2REAL8Vector(py_obj, vec, __UNIQUE_ID(_id))
#define _PyArray2REAL8Vector(obj, vec, x) \
PyArrayObject *__PASTE(x,parr) = (PyArrayObject *)PyArray_FromAny(obj, PyArray_DescrFromType(NPY_DOUBLE), 1, 0, NPY_ARRAY_DEFAULT | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST, NULL); \
if (!__PASTE(x,parr)) return NULL; \
INT __PASTE(x,ndim) = PyArray_NDIM(__PASTE(x,parr)); \
INT __PASTE(x,size) = PyArray_SIZE(__PASTE(x,parr)); \
INT __PASTE(x,nptr) = PyArray_DIM(__PASTE(x,parr), __PASTE(x,ndim) - 1); \
if (__PASTE(x,size) > __PASTE(x,nptr)) \
{ \
print_err("Warning -%s: Convert ndim = %d to REAL8Vector may cause error.\n", __func__, __PASTE(x,ndim)); \
} \
REAL8 *__PASTE(x,dptr) = (REAL8 *)PyArray_DATA(__PASTE(x,parr)); \
REAL8Vector *vec = CreateREAL8Vector(__PASTE(x,nptr)); \
memcpy(vec->data, __PASTE(x,dptr), __PASTE(x,nptr) * sizeof(REAL8)); \
Py_DECREF(__PASTE(x,parr));



#define MACRO_PyArray2CHARVector(py_obj, vec) _PyArray2CHARVector(py_obj, vec, __UNIQUE_ID(_id))
#define _PyArray2CHARVector(obj, vec, x) \
PyArrayObject *__PASTE(x,parr) = (PyArrayObject *)PyArray_FromAny(obj, PyArray_DescrFromType(NPY_STRING), 1, 0, NPY_ARRAY_DEFAULT | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST , NULL); \
if (!__PASTE(x,parr)) return NULL; \
INT __PASTE(x,ndim) = PyArray_NDIM(__PASTE(x,parr)); \
INT __PASTE(x,size) = PyArray_SIZE(__PASTE(x,parr)); \
INT __PASTE(x,nptr) = PyArray_DIM(__PASTE(x,parr), PyArray_NDIM(__PASTE(x,parr)) - 1); \
if (__PASTE(x,size) > __PASTE(x,nptr)) \
{ \
print_err("Warning -%s: Convert ndim = %d to CHARVector may cause error.\n", __func__, __PASTE(x,ndim)); \
} \
UINT __PASTE(x,i); \
npy_intp __PASTE(x,stride) = PyArray_STRIDE(__PASTE(x,parr), 0); \
const CHAR * __PASTE(x,dataptr) = PyArray_BYTES(__PASTE(x,parr)); \
PyObject *__PASTE(x,sngl_obj); \
CHARVector *vec = CreateCHARVector(__PASTE(x,nptr), __PASTE(x,stride)); \
for(__PASTE(x,i)=0; __PASTE(x,i)<__PASTE(x,nptr); ++__PASTE(x,i)) \
{ \
__PASTE(x,sngl_obj) = PyArray_GETITEM(__PASTE(x,parr), __PASTE(x,dataptr) + __PASTE(x,i) * __PASTE(x,stride)); \
memcpy(vec->data[__PASTE(x,i)], PyBytes_AsString(__PASTE(x,sngl_obj)), __PASTE(x,stride)); \
} \
Py_DECREF(__PASTE(x,parr));

#define MACRO_PyArray2REAL8Array2DIM(py_obj, arr) _PyArray2REAL8Array2DIM(py_obj, arr, __UNIQUE_ID(_id))
#define _PyArray2REAL8Array2DIM(obj, arr, x) \
PyArrayObject *__PASTE(x,data) = (PyArrayObject *)PyArray_FromAny(obj, PyArray_DescrFromType(NPY_DOUBLE), 1, 0, NPY_ARRAY_DEFAULT | NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST, NULL); \
if (!__PASTE(x,data)) return NULL; \
if (PyArray_NDIM(__PASTE(x,data)) != 2) return NULL; \
INT __PASTE(x,ns2) = PyArray_DIM(__PASTE(x,data), 1); \
INT __PASTE(x,ns1) = PyArray_DIM(__PASTE(x,data), 0); \
REAL8Array *arr = CreateREAL8Array(2, __PASTE(x,ns1), __PASTE(x,ns2)); \
REAL8 *__PASTE(x,dptr) = (REAL8 *)PyArray_DATA(__PASTE(x,data)); \
memcpy(arr->data, __PASTE(x,dptr), arr->size * sizeof(REAL8)); \
Py_DECREF(__PASTE(x,data));



#define MACRO_REAL8Array2PyArray(arr, ret) _REAL8Array2PyArray(arr, ret, __UNIQUE_ID(_id))
#define _REAL8Array2PyArray(arr, ret, x) \
INT __PASTE(x,ndim) = arr->dimLength->length; \
npy_intp __PASTE(x,dims)[__PASTE(x,ndim)]; \
UINT __PASTE(x,i); \
for (__PASTE(x,i)=0; __PASTE(x,i) < __PASTE(x,ndim); ++__PASTE(x,i)) \
{ \
__PASTE(x,dims)[__PASTE(x,i)] = arr->dimLength->data[__PASTE(x,i)]; \
} \
PyArrayObject *ret = (PyArrayObject *)PyArray_Empty(__PASTE(x,ndim), __PASTE(x,dims), PyArray_DescrFromType(NPY_DOUBLE), 0); \
REAL8 *__PASTE(x,rptr) = (REAL8 *)PyArray_DATA(ret); \
memcpy(__PASTE(x,rptr), arr->data, arr->size * sizeof(REAL8)); \
DestroyREAL8Array(arr);



#endif

