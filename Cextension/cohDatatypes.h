/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_COHDATATYPES__
#define __INCLUDE_COHDATATYPES__

#include <math.h>
#include <complex.h>
#include <stdio.h>


typedef short SINT;
typedef long LINT;
typedef int INT;
typedef unsigned UINT;

typedef char CHAR;

typedef float REAL4;
typedef double REAL8;
typedef long double REAL16;

typedef double complex COMPLEX16;


#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif

/** < STRUCTURE DEF > **/
typedef struct tagUINTVector
{
    UINT length;
    UINT *data;
}
UINTVector;


typedef struct tagREAL8Vector
{
    UINT length; /** Number of elements **/
    REAL8 *data; /** Pointer to the data array **/
}
REAL8Vector;

typedef struct tagREAL8Array
{
    UINTVector *dimLength; /** Vector of array dimensions **/
    REAL8 *data; /** Pointer to the data array **/
    UINT size; /** Number of data **/
}
REAL8Array;

/** < CONTROL ERROR VALUE DEF (CEV) > **/
enum ControlErrorValue
{
    CEV_SUCCESS     = 0, /** Success return value **/
    CEV_FAILURE     = 1, /** Failure return **/
};


#endif

