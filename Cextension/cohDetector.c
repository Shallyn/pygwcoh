/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "cohDetector.h"

const LALUnit lalStrainUnit        = {  0, { 0, 0, 0, 0, 0, 1, 0}, { 0, 0, 0, 0, 0, 0, 0} };
const LALDetector lalCachedDetectors[LAL_NUM_DETECTORS] = {
    LAL_DETECTOR_STRUCT( TAMA_300, IFODIFF ),
    LAL_DETECTOR_STRUCT( VIRGO, IFODIFF ),
    LAL_DETECTOR_STRUCT( GEO_600, IFODIFF ),
    LAL_DETECTOR_STRUCT( LHO_2K, IFODIFF ),
    LAL_DETECTOR_STRUCT( LHO_4K, IFODIFF ),
    LAL_DETECTOR_STRUCT( LLO_4K, IFODIFF ),
    LAL_DETECTOR_STRUCT( CIT_40, IFODIFF ),
    LAL_DETECTOR_STRUCT( ALLEGRO_320, CYLBAR ),
    LAL_DETECTOR_STRUCT( AURIGA, CYLBAR ),
    LAL_DETECTOR_STRUCT( EXPLORER, CYLBAR ),
    LAL_DETECTOR_STRUCT( NIOBE, CYLBAR ),
    LAL_DETECTOR_STRUCT( NAUTILUS, CYLBAR ),
    LAL_DETECTOR_STRUCT( ET1, IFODIFF ),
    LAL_DETECTOR_STRUCT( ET2, IFODIFF ),
    LAL_DETECTOR_STRUCT( ET3, IFODIFF ),
    LAL_DETECTOR_STRUCT( ET0, IFODIFF ),
    LAL_DETECTOR_STRUCT( KAGRA, IFODIFF ),
    LAL_DETECTOR_STRUCT( LIO_4K, IFODIFF ),
};


void antenna_pattern_comm(const double    D[3][3],
                          const double    ra,
                          const double    dec,
                          const double    psi,
                          const double    gmst,
                          double    *fplus,
                          double    *fcross)
{
    const double gha = gmst - ra;
    int i;
    double X[3];
    double Y[3];
    /* pre-compute trig functions */
    const double cosgha = cos(gha);
    const double singha = sin(gha);
    const double cosdec = cos(dec);
    const double sindec = sin(dec);
    const double cospsi = cos(psi);
    const double sinpsi = sin(psi);
    
    /* Eq. (B4) of [ABCF].  Note that dec = pi/2 - theta, and gha =
     * -phi where theta and phi are the standard spherical coordinates
     * used in that paper. */
    X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
    X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
    X[2] =  sinpsi * cosdec;
    
    /* Eq. (B5) of [ABCF].  Note that dec = pi/2 - theta, and gha =
     * -phi where theta and phi are the standard spherical coordinates
     * used in that paper. */
    Y[0] =  sinpsi * singha - cospsi * cosgha * sindec;
    Y[1] =  sinpsi * cosgha + cospsi * singha * sindec;
    Y[2] =  cospsi * cosdec;
    
    /* Now compute Eq. (B7) of [ABCF] for each polarization state, i.e.,
     * with s+=1 and sx=0 to get F+, with s+=0 and sx=1 to get Fx */
    *fplus = *fcross = 0.0;
    for(i = 0; i < 3; i++) {
        const double DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
        const double DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
        *fplus  += X[i] * DX - Y[i] * DY;
        *fcross += X[i] * DY + Y[i] * DX;
    }
}

double time_delay_from_earth_center(const double    detloc[3],
                                    const double          ra,
                                    const double          de,
                                    const double    gpstime)
{
    return TimeDelayFromEarthCenter(detloc, ra, de, gpstime);
}


void set_detector(const char *prefix, const LALDetector **detector)
{
    int d;
    for (d = 0; d < LAL_NUM_DETECTORS; ++d)
    {
        if (strcmp(prefix, lalCachedDetectors[d].frDetector.prefix) == 0)
        {
            *detector = lalCachedDetectors + d;
            break;
        }
    }
}

