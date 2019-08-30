/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_COHCORE__
#define __INCLUDE_COHCORE__

#include "pyUtils.h"
#include "cohDetector.h"

INT Gpc_and_delay_interpolate(CHAR *ifo,
                              REAL8Vector *ra,
                              REAL8Vector *de,
                              REAL8Vector *SNR,
                              REAL8Vector *time_SNR,
                              REAL8Vector *time,
                              REAL8 geocent_time,
                              REAL8Array **Gpc,
                              REAL8Array **SNR_matrix);


#endif

