
#include "pyUtils.h"



// Core
static PyObject * Fun_Gpc_time_pix(PyObject *self, PyObject *args)
{
    PyObject *py_SNR_real, *py_SNR_imag, *py_SNR_time;
    PyObject *py_ra, *py_de, *py_time;
    CHAR *ifo;
    REAL8 geocent_time;
    if(!PyArg_ParseTuple(args, "OOOOOOsd",
                         &py_SNR_real, &py_SNR_imag, &py_SNR_time,
                         &py_ra, &py_de, &py_time,
                         &ifo, &geocent_time))
    {
        return NULL;
    }
    // Wrapping
    MACRO_PyArray2REAL8Vector(py_ra, vec_ra);
    MACRO_PyArray2REAL8Vector(py_de, vec_de);
    MACRO_PyArray2REAL8Vector(py_SNR_real, vec_SNR_real);
    MACRO_PyArray2REAL8Vector(py_SNR_imag, vec_SNR_imag);
    MACRO_PyArray2REAL8Vector(py_SNR_time, vec_SNR_time);
    MACRO_PyArray2REAL8Vector(py_time, vec_time);
    
    REAL8Array *Garr = NULL;
    REAL8Array *Sarr_real = NULL, *Sarr_imag = NULL;
    Py_BEGIN_ALLOW_THREADS;
    NPY_SIGINT_ON;
    if (Gpc_and_delay_interpolate(ifo, vec_ra, vec_de,
                                  vec_SNR_real, vec_SNR_imag, vec_SNR_time,
                                  vec_time, geocent_time,
                                  &Garr, &Sarr_real, &Sarr_imag) != CEV_SUCCESS)
    {
        print_err("Error : Fail to calculate Gpc\n");
    }
    NPY_SIGINT_OFF;
    Py_END_ALLOW_THREADS;
    
    MACRO_REAL8Array2PyArray(Garr, Gret);
    MACRO_REAL8Array2PyArray(Sarr_real, Sret_real);
    MACRO_REAL8Array2PyArray(Sarr_imag, Sret_imag);
    DestroyREAL8Vector(vec_SNR_real);
    DestroyREAL8Vector(vec_SNR_imag);
    DestroyREAL8Vector(vec_SNR_time);
    DestroyREAL8Vector(vec_ra);
    DestroyREAL8Vector(vec_de);
    DestroyREAL8Vector(vec_time);
    
    PyObject *out;
    out = Py_BuildValue("OOO", Gret, Sret_real, Sret_imag);
    return out;
}




static PyMethodDef PyGWCOH_Methods[] = {
    {"Gpc_time_pix", Fun_Gpc_time_pix, METH_VARARGS},
    {NULL, NULL},
};

static struct PyModuleDef PyGWCOHModule = {
    PyModuleDef_HEAD_INIT,
    "PyGWCOH",
    NULL,
    -1,
    PyGWCOH_Methods
};

PyMODINIT_FUNC PyInit_PyGWCOH(void)
{
    PyObject *m = PyModule_Create(&PyGWCOHModule);
    import_array();
    return m;
}
