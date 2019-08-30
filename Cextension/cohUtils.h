/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_COHUTILS__
#define __INCLUDE_COHUTILS__

#include "cohDatatypes.h"
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/** < MACRO DEF > **/
/* Internal */
#define ___PASTE(a,b) a##b
#define __PASTE(a,b) ___PASTE(a,b)
#define __UNIQUE_ID(prefix) __PASTE(__PASTE(__UNIQUE_ID_, prefix), __COUNTER__)

#define _SWAP(a,b,tmp) \
do{typeof(a) tmp = (a);(a)=(b);(b)=tmp;} while(0)

#define SWAP(a,b) \
_SWAP(a,b,__UNIQUE_ID(tmp))

/** < FUNCTIONS DECLARATION > **/

/* Used for standart print */
INT print_err(const CHAR *fmt, ...);

/* Used for creating & destroying basic data struct */
UINTVector* CreateUINTVector(UINT length);
void DestroyUINTVector(UINTVector* vector);

REAL8Vector* CreateREAL8Vector(UINT length);
void DestroyREAL8Vector(REAL8Vector* vector);

REAL8Array *CreateREAL8Array_comm(UINT *dims, UINT ndim);
REAL8Array *CreateREAL8Array(UINT ndim, ...);
void DestroyREAL8Array(REAL8Array* array);

/* REAL8Array */
INT getArrayValue_comm(REAL8 *val, REAL8Array *arr, UINT *cood);
INT getArrayValue(REAL8 *val, REAL8Array *arr, UINT ndim,...);
void setArrayValue_comm(REAL8 val, REAL8Array *arr, UINT *cood);
INT setArrayValue(REAL8 val, REAL8Array *arr, UINT ndim,...);


#endif

