/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "cohCore.h"
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

INT Gpc_and_delay_interpolate(CHAR *ifo,
                              REAL8Vector *ra,
                              REAL8Vector *de,
                              REAL8Vector *SNR_real,
                              REAL8Vector *SNR_imag,
                              REAL8Vector *time_SNR,
                              REAL8Vector *time,
                              REAL8 geocent_time,
                              REAL8Array **Gpc_matrix,
                              REAL8Array **SNR_real_matrix,
                              REAL8Array **SNR_imag_matrix)
{
    if (ra->length != de->length)
    {
        print_err("Error -%s: Length of ra and de is different.\n",__func__);
        return CEV_FAILURE;
    }
    if( (SNR_real->length != time_SNR->length) || (SNR_imag->length != time_SNR->length))
    {
        print_err("Error -%s: Length of sigma2 and ifo is different.\n", __func__);
        return CEV_FAILURE;
    }
    UINT npix = ra->length;
    UINT ntime = time->length;
    REAL8Array *arr_Gpc = CreateREAL8Array(2, npix, 2);
    REAL8Array *arr_SNR_real = CreateREAL8Array(2, ntime, npix);
    REAL8Array *arr_SNR_imag = CreateREAL8Array(2, ntime, npix);

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline_real = gsl_spline_alloc (gsl_interp_cspline, SNR_real->length);
    gsl_spline *spline_imag = gsl_spline_alloc (gsl_interp_cspline, SNR_imag->length);
    gsl_spline_init (spline_real, time_SNR->data, SNR_real->data, SNR_real->length);
    gsl_spline_init (spline_imag, time_SNR->data, SNR_imag->data, SNR_imag->length);

    
    REAL8 ra_i, de_i;
    REAL8 fplus, fcross;
    REAL8 gmst, delay;
    const LALDetector *detector;
    set_detector(ifo, &detector);
    UINT i,j, cood[2];
    REAL8 snr_itp;
    for(i=0; i < npix; ++i)
    {
        cood[0] = i;
        ra_i = ra->data[i];
        de_i = de->data[i];
        gmst = GreenwichMeanSiderealTime(geocent_time);
        antenna_pattern_comm(detector->response, ra_i, de_i, 0.0, gmst, &fplus, &fcross);
        delay = time_delay_from_earth_center(detector->location, ra_i, de_i, geocent_time);
        cood[1] = 0;
        setArrayValue_comm(fplus, arr_Gpc, cood);
        cood[1] = 1;
        setArrayValue_comm(fcross, arr_Gpc, cood);
        for(j = 0; j < ntime; ++j)
        {
            cood[0] = j;
            cood[1] = i;
            gsl_spline_eval_e (spline_real, time->data[j] + delay, acc, &snr_itp);
            setArrayValue_comm(snr_itp, arr_SNR_real, cood);
            gsl_spline_eval_e (spline_imag, time->data[j] + delay, acc, &snr_itp);
            setArrayValue_comm(snr_itp, arr_SNR_imag, cood);
        }
    }
    gsl_spline_free (spline_real);
    gsl_spline_free (spline_imag);
    gsl_interp_accel_free (acc);

    *Gpc_matrix = arr_Gpc;
    *SNR_real_matrix = arr_SNR_real;
    *SNR_imag_matrix = arr_SNR_imag;
    return CEV_SUCCESS;
}


