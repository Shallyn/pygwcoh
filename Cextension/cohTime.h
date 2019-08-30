/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_COHTIME__
#define __INCLUDE_COHTIME__
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

typedef int64_t gwINT8;

// Constants & Macro
#define LAL_PI        3.1415926535897932384626433832795029  /**< pi */
#define LAL_C_SI      299792458 /**< Speed of light in vacuo, m s^-1 */

#define BILLION_INT8 INT64_C( 1000000000 )

#define EPOCH_UNIX_GPS 315964800

#define EPOCH_J2000_0_JD 2451545.0         /**< Julian Day (UTC) of the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define EPOCH_J2000_0_TAI_UTC 32           /**< Leap seconds (TAI-UTC) on the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define EPOCH_J2000_0_GPS 630763213        /**< GPS seconds of the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define EPOCH_GPS_JD 2444244.5             /**< Julian Day (UTC) of the GPS epoch (1980 JAN 6 0h UTC) */
#define EPOCH_GPS_TAI_UTC 19               /**< Leap seconds (TAI-UTC) on the GPS epoch (1980 JAN 6 0h UTC) */
#define MJD_REF 2400000.5                  /**< Reference Julian Day for Mean Julian Day. */
#define JD_TO_MJD(jd) ((jd) - MJD_REF) /**< Modified Julian Day for specified civil time structure. */


// Datatypes
typedef struct
tagGPSTime
{
    int gpsSeconds; /**< Seconds since 0h UTC 6 Jan 1980. */
    int gpsNanoSeconds; /**< Residual nanoseconds. */
}
GPSTime;

enum {
    SUCCESS = 0,
    FAILURE = -1
};

/* Functions */
double GreenwichMeanSiderealTime(const double gpstime);
double TimeDelayFromEarthCenter(const double detector_earthfixed_xyz_metres[3], double source_right_ascension_radians, double source_declination_radians, const double gpstime);
double GPStoREAL8(GPSTime *gps);

#endif

