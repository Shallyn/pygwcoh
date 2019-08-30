/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "cohTime.h"

static const struct leaps_table { double jd; int gpssec; int taiutc; } leaps[] =
{
    {2444239.5,    -43200, 19},  /* 1980-Jan-01 */
    {2444786.5,  46828800, 20},  /* 1981-Jul-01 */
    {2445151.5,  78364801, 21},  /* 1982-Jul-01 */
    {2445516.5, 109900802, 22},  /* 1983-Jul-01 */
    {2446247.5, 173059203, 23},  /* 1985-Jul-01 */
#if 0
    /* NOTE: IF THIS WERE A NEGATIVE LEAP SECOND, INSERT AS FOLLOWS */
    {2447161.5, 252028803, 22},  /* 1988-Jan-01 EXAMPLE ONLY! */
#endif
    {2447161.5, 252028804, 24},  /* 1988-Jan-01 */
    {2447892.5, 315187205, 25},  /* 1990-Jan-01 */
    {2448257.5, 346723206, 26},  /* 1991-Jan-01 */
    {2448804.5, 393984007, 27},  /* 1992-Jul-01 */
    {2449169.5, 425520008, 28},  /* 1993-Jul-01 */
    {2449534.5, 457056009, 29},  /* 1994-Jul-01 */
    {2450083.5, 504489610, 30},  /* 1996-Jan-01 */
    {2450630.5, 551750411, 31},  /* 1997-Jul-01 */
    {2451179.5, 599184012, 32},  /* 1999-Jan-01 */
    {2453736.5, 820108813, 33},  /* 2006-Jan-01 */
    {2454832.5, 914803214, 34},  /* 2009-Jan-01 */
    {2456109.5, 1025136015, 35}, /* 2012-Jul-01 */
    {2457204.5, 1119744016, 36}, /* 2015-Jul-01 */
    {2457754.5, 1167264017, 37}, /* 2017-Jan-01 */
};
static const int numleaps = sizeof( leaps ) / sizeof( *leaps );

static double dotprod(const double vec1[3], const double vec2[3])
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}


/* GPSTime handler */
static int delta_tai_utc( int gpssec )
{
    int leap;
    
    /* assume calling function has already checked this */
    /*
     if ( gpssec <= leaps[0].gpssec )
     {
     fprintf( stderr, "error: don't know leap seconds before gps time %d\n",
     leaps[0].gpssec );
     abort();
     }
     */
    
    for ( leap = 1; leap < numleaps; ++leap )
        if ( gpssec == leaps[leap].gpssec )
            return leaps[leap].taiutc - leaps[leap-1].taiutc;
    
    return 0;
}








/* GreenWitchSIderealTime */
double ConvertCivilTimeToJD( const struct tm *civil)
{
    const int sec_per_day = 60 * 60 * 24; /* seconds in a day */
    int year, month, day, sec;
    double jd;
    
    /* this routine only works for dates after 1900 */
    if ( civil->tm_year <= 0 )
    {
        fprintf(stderr, "XLAL Error - Year must be after 1900\n" );
        return (double)FAILURE;
    }
    
    year  = civil->tm_year + 1900;
    month = civil->tm_mon + 1;     /* month is in range 1-12 */
    day   = civil->tm_mday;        /* day is in range 1-31 */
    sec   = civil->tm_sec + 60*(civil->tm_min + 60*(civil->tm_hour)); /* seconds since midnight */
    
    jd = 367*year - 7*(year + (month + 9)/12)/4 + 275*month/9 + day + 1721014;
    /* note: Julian days start at noon: subtract half a day */
    jd += (double)sec/(double)sec_per_day - 0.5;
    
    return jd;
} // XLALConvertCivilTimeToJD()



int LeapSeconds( int gpssec )
{
    int leap;
    
    if ( gpssec < leaps[0].gpssec )
    {
        fprintf(stderr, "XLAL Error - Don't know leap seconds before GPS time %d\n",
                leaps[0].gpssec );
        return -1;
    }
    
    /* scan leap second table and locate the appropriate interval */
    for ( leap = 1; leap < numleaps; ++leap )
        if ( gpssec < leaps[leap].gpssec )
            break;
    
    return leaps[leap-1].taiutc;
}



struct tm * GPSToUTC(struct tm *utc, int gpssec)
{
    time_t unixsec;
    int leapsec;
    int delta;
    leapsec = LeapSeconds( gpssec );
    if ( leapsec < 0 )
        return NULL;
    unixsec  = gpssec - leapsec + EPOCH_GPS_TAI_UTC; /* get rid of leap seconds */
    unixsec += EPOCH_UNIX_GPS; /* change to unix epoch */
    memset( utc, 0, sizeof( *utc ) ); /* blank out utc structure */
    utc = gmtime_r( &unixsec, utc );
    /* now check to see if we need to add a 60th second to UTC */
    if ( ( delta = delta_tai_utc( gpssec ) ) > 0 )
        utc->tm_sec += 1; /* delta only ever is one, right?? */
    return utc;
}


double GreenwichSiderealTime(const double gpstime, double equation_of_equinoxes)
{
    struct tm utc;
    double julian_day;
    double t_hi, t_lo;
    double t;
    double sidereal_time;
    
    /*
     * Convert GPS seconds to UTC.  This is where we pick up knowledge
     * of leap seconds which are required for the mapping of atomic
     * time scales to celestial time scales.  We deal only with integer
     * seconds.
     */
    double gpssec;
    double gpsnanosec = modf(gpstime, &gpssec);
    
    if(!GPSToUTC(&utc, gpssec))
        return (double)FAILURE;
    
    /*
     * And now to Julian day number.  Again, only accurate to integer
     * seconds.
     */
    
    julian_day = ConvertCivilTimeToJD(&utc);
    
    /*
     * Convert Julian day number to the number of centuries since the
     * Julian epoch (1 century = 36525.0 days).  Here, we incorporate
     * the fractional part of the seconds.  For precision, we keep
     * track of the most significant and least significant parts of the
     * time separately.  The original code in NOVAS-C determined t_hi
     * and t_lo from Julian days, with t_hi receiving the integer part
     * and t_lo the fractional part.  Because LAL's Julian day routine
     * is accurate to the second, here the hi/lo split is most
     * naturally done at the integer seconds boundary.  Note that the
     * "hi" and "lo" components have the same units and so the split
     * can be done anywhere.
     */
    
    t_hi = (julian_day - EPOCH_J2000_0_JD) / 36525.0;
    
    t_lo = (gpsnanosec * BILLION_INT8) / (1e9 * 36525.0 * 86400.0);
    
    /*
     * Compute sidereal time in sidereal seconds.  (magic)
     */
    
    t = t_hi + t_lo;
    
    sidereal_time = equation_of_equinoxes + (-6.2e-6 * t + 0.093104) * t * t + 67310.54841;
    sidereal_time += 8640184.812866 * t_lo;
    sidereal_time += 3155760000.0 * t_lo;
    sidereal_time += 8640184.812866 * t_hi;
    sidereal_time += 3155760000.0 * t_hi;
    
    /*
     * Return radians (2 pi radians in 1 sidereal day = 86400 sidereal
     * seconds).
     */
    
    return sidereal_time * LAL_PI / 43200.0;
}



double GreenwichMeanSiderealTime(const double gpstime)
{
    return GreenwichSiderealTime(gpstime, 0.0);
}

/* Time Delay */
double ArrivalTimeDiff(
                       const double detector1_earthfixed_xyz_metres[3],
                       const double detector2_earthfixed_xyz_metres[3],
                       const double source_right_ascension_radians,
                       const double source_declination_radians,
                       const double gpstime
                       )
{
    double delta_xyz[3];
    double ehat_src[3];
    const double greenwich_hour_angle = GreenwichMeanSiderealTime(gpstime) - source_right_ascension_radians;
    
    
    /*
     * compute the unit vector pointing from the geocenter to the
     * source
     */
    
    ehat_src[0] = cos(source_declination_radians) * cos(greenwich_hour_angle);
    ehat_src[1] = cos(source_declination_radians) * -sin(greenwich_hour_angle);
    ehat_src[2] = sin(source_declination_radians);
    
    /*
     * position of detector 2 with respect to detector 1
     */
    
    delta_xyz[0] = detector2_earthfixed_xyz_metres[0] - detector1_earthfixed_xyz_metres[0];
    delta_xyz[1] = detector2_earthfixed_xyz_metres[1] - detector1_earthfixed_xyz_metres[1];
    delta_xyz[2] = detector2_earthfixed_xyz_metres[2] - detector1_earthfixed_xyz_metres[2];
    
    /*
     * Arrival time at detector 1 - arrival time at detector 2.  This
     * is positive when the wavefront arrives at detector 1 after
     * detector 2 (and so t at detector 1 is greater than t at detector
     * 2).
     */
    
    return dotprod(ehat_src, delta_xyz) / LAL_C_SI;
}


double TimeDelayFromEarthCenter(const double detector_earthfixed_xyz_metres[3], const double source_right_ascension_radians, const double source_declination_radians, const double gpstime)
{
    static const double earth_center[3] = {0.0, 0.0, 0.0};
    
    /*
     * This is positive when the wavefront arrives at the detector
     * after arriving at the geocentre.
     */
    
    return ArrivalTimeDiff(detector_earthfixed_xyz_metres, earth_center, source_right_ascension_radians, source_declination_radians, gpstime);
}

double GPStoREAL8(GPSTime *gps)
{
    return (double)gps->gpsSeconds + (double)gps->gpsNanoSeconds / (double)BILLION_INT8;
}


