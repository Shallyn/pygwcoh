/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "cohUtils.h"

/** < FUNCTIONS > **/

/* Used to print error standard error message and output message */
static INT myPrintVstderr(const CHAR *fmt, va_list ap)
{
    return vfprintf(stderr, fmt, ap);
}


INT print_err(const CHAR *fmt, ...)
{
    INT n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = myPrintVstderr(fmt, ap);
    va_end(ap);
    return n;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

/** < create & destroy struct > **/
/* UINTVector */
UINTVector* CreateUINTVector(UINT length)
{
    UINTVector* vector;
    vector = (UINTVector*)malloc(sizeof(*vector));
    vector->length = length;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (UINT *)malloc(length * sizeof( *vector->data));
    }
    
    return vector;
}
/* Destroy */
void DestroyUINTVector(UINTVector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}

/* REAL8Vector */
REAL8Vector* CreateREAL8Vector(UINT length)
{
    REAL8Vector* vector;
    vector = (REAL8Vector*)malloc(sizeof(*vector));
    vector->length = length;
    if ( ! length ) /* zero length: set data pointer to be NULL */
    {
        vector->data = NULL;
    }
    else /* non-zero length: allocate memory for data */
    {
        vector->data = (REAL8 *)malloc(length * sizeof( *vector->data));
    }
    
    return vector;
}

/* Destroy */
void DestroyREAL8Vector(REAL8Vector* vector)
{
    if(NULL == vector)
    {
        return;
    }
    if(vector->data)
    {
        vector->length = 0;
        free(vector->data);
    }
    vector->data = NULL;
    free(vector);
    vector = NULL;
    return;
}

/* REAL8Array */
REAL8Array *CreateREAL8Array_comm(UINT *dims, UINT ndim)
{
    UINTVector dimLength;
    UINT size = 1;
    REAL8Array *arr;
    dimLength.length = ndim;
    dimLength.data = dims;
    
    UINT i;
    for(i = 0; i < ndim; ++i)
        size *= dimLength.data[i];
    arr = (REAL8Array *)malloc(sizeof(*arr));
    arr->size = size;
    arr->dimLength = CreateUINTVector(ndim);
    if (! arr->dimLength)
    {
        free(arr);
        arr = NULL;
        return NULL;
    }
    
    memcpy(arr->dimLength->data, dimLength.data, ndim*sizeof(*arr->dimLength->data));
    
    arr->data = (REAL8 *)malloc( size * sizeof(*arr->data));
    if(!arr->data)
    {
        DestroyUINTVector(arr->dimLength);
        free(arr);
        arr = NULL;
        return NULL;
    }
    return arr;
    
}

REAL8Array *CreateREAL8Array(UINT ndim, ...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return NULL;
    }
    va_list ap;
    UINT dims[ndim];
    UINT dim;
    
    va_start(ap, ndim);
    for(dim = 0; dim<ndim; ++dim)
        dims[dim] = va_arg(ap, UINT);
    va_end(ap);
    return CreateREAL8Array_comm(dims, ndim);
}

/* Destroy */
void DestroyREAL8Array(REAL8Array* array)
{
    DestroyUINTVector( array->dimLength );
    free( array->data );
    array->data = NULL;
    free( array );
    array = NULL;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */


/** < REAL8Array > **/
/* Array element: for [i,j,k],
 * get corresponding dim from Array->dimLength [0],[1],[2],
 * data = Array->data[i + j*di + k*di*dj]
 */

static inline UINT Cood2Index(UINT *cood, UINT *shape, UINT ndim)
{
    UINT index = 0, i;
    UINT cum = 1;
    for(i=0; i < ndim; ++i)
    {
        cum *= i == 0 ? 1 : shape[ndim - i];
        index += cood[ndim - i - 1] * cum;
    }
    return index;
}

static inline INT Index2Cood(UINT **cood ,UINT index, UINT *shape, UINT ndim, UINT size)
{
    // The length of cood should be ndim, be careful!!
    UINT cum = size, i, tmp;
    for(i=0; i < ndim - 1; ++i)
    {
        cum /= shape[i];
        tmp = index / cum;
        (*cood)[i] = tmp;
        index -= cum * tmp;
    }
    (*cood)[ndim - 1] = index;
    return CEV_SUCCESS;
}

/* get value for REAL8Array */
INT getArrayValue_comm(REAL8 *val, REAL8Array *arr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, arr->dimLength->data, arr->dimLength->length);
    *val = arr->data[index];
    return CEV_SUCCESS;
}

INT getArrayValue(REAL8 *val, REAL8Array *arr, UINT ndim,...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return CEV_FAILURE;
    }
    
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != arr->dimLength->length)
    {
        print_err("Error %s: The shape of input cood is not equal to the shape of the input array.\n", __func__);
        return CEV_FAILURE;
    }
    return getArrayValue_comm(val, arr, dims);
}

/* set value for REAL8Array */
void setArrayValue_comm(REAL8 val, REAL8Array *arr, UINT *cood)
{
    UINT index;
    index = Cood2Index(cood, arr->dimLength->data, arr->dimLength->length);
    arr->data[index] = val;
}

INT setArrayValue(REAL8 val, REAL8Array *arr, UINT ndim,...)
{
    enum { maxdim = 16 };
    if(ndim > maxdim)
    {
        print_err("Error %s: dimension input %d is greater than max dimension %d", __func__, ndim, maxdim);
        return CEV_FAILURE;
    }
    UINT dims[ndim];
    UINT i;
    
    va_list ap;
    va_start(ap, ndim);
    for(i=0; i < ndim; ++i)
    {
        dims[i] = va_arg(ap, UINT);
    }
    va_end(ap);
    
    if (ndim != arr->dimLength->length)
    {
        print_err("Error %s: The shape of input cood is not equal to the shape of the input array.\n", __func__);
        return CEV_FAILURE;
    }
    
    setArrayValue_comm(val, arr, dims);
    return CEV_SUCCESS;
}

