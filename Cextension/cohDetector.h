/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#ifndef __INCLUDE_COHDETECTOR__
#define __INCLUDE_COHDETECTOR__

#include "cohTime.h"

enum {
    LALUnitIndexMeter, /**< The meter index. */
    LALUnitIndexKiloGram, /**< The kilogram index. */
    LALUnitIndexSecond, /**< The second index. */
    LALUnitIndexAmpere, /**< The ampere index. */
    LALUnitIndexKelvin, /**< The kelvin index. */
    LALUnitIndexStrain, /**< The strain index. */
    LALUnitIndexADCCount, /**< The ADC counts index. */
    LALNumUnits         /**< The number of units. */
};

/** Enumeration of Detectors: follows order of DQ bit assignments */
enum {
    LAL_TAMA_300_DETECTOR    =    0,
    LAL_VIRGO_DETECTOR    =    1,
    LAL_GEO_600_DETECTOR    =    2,
    LAL_LHO_2K_DETECTOR    =    3,
    LAL_LHO_4K_DETECTOR    =    4,
    LAL_LLO_4K_DETECTOR    =    5,
    LAL_CIT_40_DETECTOR    =    6,
    LAL_ALLEGRO_DETECTOR    =    7,
    LAL_AURIGA_DETECTOR    =    8,
    LAL_EXPLORER_DETECTOR    =    9,
    LAL_NIOBE_DETECTOR    =    10,
    LAL_NAUTILUS_DETECTOR    =    11,
    LAL_ET1_DETECTOR    =    12,
    LAL_ET2_DETECTOR    =    13,
    LAL_ET3_DETECTOR    =    14,
    LAL_ET0_DETECTOR    =    15,
    LAL_KAGRA_DETECTOR    =    16,
    LAL_LIO_4K_DETECTOR =   17,
    LAL_NUM_DETECTORS    =    18
};


typedef struct tagLALUnit {
    int16_t powerOfTen; /**< Overall power-of-ten scaling is 10^\c powerOfTen. */
    int16_t unitNumerator[LALNumUnits]; /**< Array of unit power numerators. */
    uint16_t unitDenominatorMinusOne[LALNumUnits]; /**< Array of unit power denominators-minus-one. */
} LALUnit;

enum enumLALNameLength { LALNameLength = 64 };

/** Detector Constants **/

#define LAL_TAMA_300_DETECTOR_NAME                   "TAMA_300"    /**< TAMA_300 detector name string */
#define LAL_TAMA_300_DETECTOR_PREFIX                 "T1"    /**< TAMA_300 detector prefix string */
#define LAL_TAMA_300_DETECTOR_LONGITUDE_RAD          2.43536359469    /**< TAMA_300 vertex longitude (rad) */
#define LAL_TAMA_300_DETECTOR_LATITUDE_RAD           0.62267336022    /**< TAMA_300 vertex latitude (rad) */
#define LAL_TAMA_300_DETECTOR_ELEVATION_SI           90    /**< TAMA_300 vertex elevation (m) */
#define LAL_TAMA_300_DETECTOR_ARM_X_AZIMUTH_RAD      4.71238898038    /**< TAMA_300 x arm azimuth (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_AZIMUTH_RAD      3.14159265359    /**< TAMA_300 y arm azimuth (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< TAMA_300 x arm altitude (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< TAMA_300 y arm altitude (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_X_MIDPOINT_SI      150.00000000000    /**< TAMA_300 x arm midpoint (m) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_MIDPOINT_SI      150.00000000000    /**< TAMA_300 y arm midpoint (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_X_SI            -3.94640899111e+06    /**< TAMA_300 x-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_Y_SI            3.36625902802e+06    /**< TAMA_300 y-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_Z_SI            3.69915069233e+06    /**< TAMA_300 z-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_ARM_X_DIRECTION_X               0.64896940530    /**< TAMA_300 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_X_DIRECTION_Y               0.76081450498    /**< TAMA_300 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_X_DIRECTION_Z               -0.00000000000    /**< TAMA_300 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_X               -0.44371376921    /**< TAMA_300 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_Y               0.37848471479    /**< TAMA_300 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_Z               -0.81232223390    /**< TAMA_300 z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name VIRGO 3km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * VIRGO 3km Interferometric Detector.
 */
/*@{*/
#define LAL_VIRGO_DETECTOR_NAME                   "VIRGO"    /**< VIRGO detector name string */
#define LAL_VIRGO_DETECTOR_PREFIX                 "V1"    /**< VIRGO detector prefix string */
#define LAL_VIRGO_DETECTOR_LONGITUDE_RAD          0.18333805213    /**< VIRGO vertex longitude (rad) */
#define LAL_VIRGO_DETECTOR_LATITUDE_RAD           0.76151183984    /**< VIRGO vertex latitude (rad) */
#define LAL_VIRGO_DETECTOR_ELEVATION_SI           51.884    /**< VIRGO vertex elevation (m) */
#define LAL_VIRGO_DETECTOR_ARM_X_AZIMUTH_RAD      0.33916285222    /**< VIRGO x arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_AZIMUTH_RAD      5.05155183261    /**< VIRGO y arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< VIRGO x arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< VIRGO y arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_MIDPOINT_SI      1500.00000000000    /**< VIRGO x arm midpoint (m) */
#define LAL_VIRGO_DETECTOR_ARM_Y_MIDPOINT_SI      1500.00000000000    /**< VIRGO y arm midpoint (m) */
#define LAL_VIRGO_VERTEX_LOCATION_X_SI            4.54637409900e+06    /**< VIRGO x-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Y_SI            8.42989697626e+05    /**< VIRGO y-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Z_SI            4.37857696241e+06    /**< VIRGO z-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_ARM_X_DIRECTION_X               -0.70045821479    /**< VIRGO x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Y               0.20848948619    /**< VIRGO y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Z               0.68256166277    /**< VIRGO z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_X               -0.05379255368    /**< VIRGO x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Y               -0.96908180549    /**< VIRGO y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Z               0.24080451708    /**< VIRGO z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name GEO 600m Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * GEO 600m Interferometric Detector.
 */
/*@{*/
#define LAL_GEO_600_DETECTOR_NAME                   "GEO_600"    /**< GEO_600 detector name string */
#define LAL_GEO_600_DETECTOR_PREFIX                 "G1"    /**< GEO_600 detector prefix string */
#define LAL_GEO_600_DETECTOR_LONGITUDE_RAD          0.17116780435    /**< GEO_600 vertex longitude (rad) */
#define LAL_GEO_600_DETECTOR_LATITUDE_RAD           0.91184982752    /**< GEO_600 vertex latitude (rad) */
#define LAL_GEO_600_DETECTOR_ELEVATION_SI           114.425    /**< GEO_600 vertex elevation (m) */
#define LAL_GEO_600_DETECTOR_ARM_X_AZIMUTH_RAD      1.19360100484    /**< GEO_600 x arm azimuth (rad) */
#define LAL_GEO_600_DETECTOR_ARM_Y_AZIMUTH_RAD      5.83039279401    /**< GEO_600 y arm azimuth (rad) */
#define LAL_GEO_600_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< GEO_600 x arm altitude (rad) */
#define LAL_GEO_600_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< GEO_600 y arm altitude (rad) */
#define LAL_GEO_600_DETECTOR_ARM_X_MIDPOINT_SI      300.00000000000    /**< GEO_600 x arm midpoint (m) */
#define LAL_GEO_600_DETECTOR_ARM_Y_MIDPOINT_SI      300.00000000000    /**< GEO_600 y arm midpoint (m) */
#define LAL_GEO_600_VERTEX_LOCATION_X_SI            3.85630994926e+06    /**< GEO_600 x-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_VERTEX_LOCATION_Y_SI            6.66598956317e+05    /**< GEO_600 y-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_VERTEX_LOCATION_Z_SI            5.01964141725e+06    /**< GEO_600 z-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_ARM_X_DIRECTION_X               -0.44530676905    /**< GEO_600 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_X_DIRECTION_Y               0.86651354130    /**< GEO_600 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_X_DIRECTION_Z               0.22551311312    /**< GEO_600 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_X               -0.62605756776    /**< GEO_600 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_Y               -0.55218609524    /**< GEO_600 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_Z               0.55058372486    /**< GEO_600 z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name LIGO Hanford Observatory 2km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 2km Interferometric Detector.
 */
/*@{*/
#define LAL_LHO_2K_DETECTOR_NAME                   "LHO_2k"    /**< LHO_2k detector name string */
#define LAL_LHO_2K_DETECTOR_PREFIX                 "H2"    /**< LHO_2k detector prefix string */
#define LAL_LHO_2K_DETECTOR_LONGITUDE_RAD          -2.08405676917    /**< LHO_2k vertex longitude (rad) */
#define LAL_LHO_2K_DETECTOR_LATITUDE_RAD           0.81079526383    /**< LHO_2k vertex latitude (rad) */
#define LAL_LHO_2K_DETECTOR_ELEVATION_SI           142.554    /**< LHO_2k vertex elevation (m) */
#define LAL_LHO_2K_DETECTOR_ARM_X_AZIMUTH_RAD      5.65487724844    /**< LHO_2k x arm azimuth (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_AZIMUTH_RAD      4.08408092164    /**< LHO_2k y arm azimuth (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_X_ALTITUDE_RAD     -0.00061950000    /**< LHO_2k x arm altitude (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00001250000    /**< LHO_2k y arm altitude (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_X_MIDPOINT_SI      1004.50000000000    /**< LHO_2k x arm midpoint (m) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_MIDPOINT_SI      1004.50000000000    /**< LHO_2k y arm midpoint (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_X_SI            -2.16141492636e+06    /**< LHO_2k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_Y_SI            -3.83469517889e+06    /**< LHO_2k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_Z_SI            4.60035022664e+06    /**< LHO_2k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_ARM_X_DIRECTION_X               -0.22389266154    /**< LHO_2k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_X_DIRECTION_Y               0.79983062746    /**< LHO_2k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_X_DIRECTION_Z               0.55690487831    /**< LHO_2k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_X               -0.91397818574    /**< LHO_2k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_Y               0.02609403989    /**< LHO_2k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_Z               -0.40492342125    /**< LHO_2k z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name LIGO Hanford Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 4km Interferometric Detector.
 */
/*@{*/
#define LAL_LHO_4K_DETECTOR_NAME                   "LHO_4k"    /**< LHO_4k detector name string */
#define LAL_LHO_4K_DETECTOR_PREFIX                 "H1"    /**< LHO_4k detector prefix string */
#define LAL_LHO_4K_DETECTOR_LONGITUDE_RAD          -2.08405676917    /**< LHO_4k vertex longitude (rad) */
#define LAL_LHO_4K_DETECTOR_LATITUDE_RAD           0.81079526383    /**< LHO_4k vertex latitude (rad) */
#define LAL_LHO_4K_DETECTOR_ELEVATION_SI           142.554    /**< LHO_4k vertex elevation (m) */
#define LAL_LHO_4K_DETECTOR_ARM_X_AZIMUTH_RAD      5.65487724844    /**< LHO_4k x arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD      4.08408092164    /**< LHO_4k y arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_ALTITUDE_RAD     -0.00061950000    /**< LHO_4k x arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00001250000    /**< LHO_4k y arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_MIDPOINT_SI      1997.54200000000    /**< LHO_4k x arm midpoint (m) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_MIDPOINT_SI      1997.52200000000    /**< LHO_4k y arm midpoint (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_X_SI            -2.16141492636e+06    /**< LHO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Y_SI            -3.83469517889e+06    /**< LHO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Z_SI            4.60035022664e+06    /**< LHO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_ARM_X_DIRECTION_X               -0.22389266154    /**< LHO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Y               0.79983062746    /**< LHO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Z               0.55690487831    /**< LHO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_X               -0.91397818574    /**< LHO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Y               0.02609403989    /**< LHO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Z               -0.40492342125    /**< LHO_4k z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name LIGO Livingston Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Livingston Observatory 4km Interferometric Detector.
 */
/*@{*/
#define LAL_LLO_4K_DETECTOR_NAME                   "LLO_4k"    /**< LLO_4k detector name string */
#define LAL_LLO_4K_DETECTOR_PREFIX                 "L1"    /**< LLO_4k detector prefix string */
#define LAL_LLO_4K_DETECTOR_LONGITUDE_RAD          -1.58430937078    /**< LLO_4k vertex longitude (rad) */
#define LAL_LLO_4K_DETECTOR_LATITUDE_RAD           0.53342313506    /**< LLO_4k vertex latitude (rad) */
#define LAL_LLO_4K_DETECTOR_ELEVATION_SI           -6.574    /**< LLO_4k vertex elevation (m) */
#define LAL_LLO_4K_DETECTOR_ARM_X_AZIMUTH_RAD      4.40317772346    /**< LLO_4k x arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD      2.83238139666    /**< LLO_4k y arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_ALTITUDE_RAD     -0.00031210000    /**< LLO_4k x arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD     -0.00061070000    /**< LLO_4k y arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_MIDPOINT_SI      1997.57500000000    /**< LLO_4k x arm midpoint (m) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_MIDPOINT_SI      1997.57500000000    /**< LLO_4k y arm midpoint (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_X_SI            -7.42760447238e+04    /**< LLO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Y_SI            -5.49628371971e+06    /**< LLO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Z_SI            3.22425701744e+06    /**< LLO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_ARM_X_DIRECTION_X               -0.95457412153    /**< LLO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Y               -0.14158077340    /**< LLO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Z               -0.26218911324    /**< LLO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_X               0.29774156894    /**< LLO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Y               -0.48791033647    /**< LLO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Z               -0.82054461286    /**< LLO_4k z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/



/**
 * \name LIGO India 4km Interferometric Detector constants
 * @warning These numbers are subject to change.
 * The following constants describe hypothetical location and geometry
 * of the LIGO India 4km Interferometric Detector that have been used
 * in several studies with LALInference. Note that these data do not
 * represent an actual prospective site.
 */
#define LAL_LIO_4K_DETECTOR_NAME                 "LIO_4k" /**< LIO_4K detector name string */
#define LAL_LIO_4K_DETECTOR_PREFIX               "I1"    /**< LIO_4K detector prefix string */
#define LAL_LIO_4K_DETECTOR_LONGITUDE_RAD        1.3340133249409993   /**< LIO_4K vertex longitude (rad; equal to 76°26') */
#define LAL_LIO_4K_DETECTOR_LATITUDE_RAD         0.2484185302005262   /**< LIO_4K vertex latitude (rad; equal to 14°14') */
#define LAL_LIO_4K_DETECTOR_ELEVATION_SI         0.0  /**< LIO_4K vertex elevation (m) */
#define LAL_LIO_4K_DETECTOR_ARM_X_AZIMUTH_RAD    1.5707963705062866   /**< LIO_4K x arm azimuth (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD    0.0   /**< LIO_4K y arm azimuth (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_X_ALTITUDE_RAD   0.0   /**< LIO_4K x arm altitude (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD   0.0   /**< LIO_4K y arm altitude (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_X_MIDPOINT_SI    2000.00000000000        /**< LIO_4K x arm midpoint (m) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_MIDPOINT_SI    2000.00000000000        /**< LIO_4K y arm midpoint (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_X_SI          1450526.82294155       /**< LIO_4K x-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_Y_SI          6011058.39047265       /**< LIO_4K y-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_Z_SI          1558018.27884102       /**< LIO_4K z-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_ARM_X_DIRECTION_X            -9.72097635269165039e-01  /**< LIO_4K x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_X_DIRECTION_Y             2.34576612710952759e-01   /**< LIO_4K y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_X_DIRECTION_Z            -4.23695567519644101e-08 /**< LIO_4K z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_X             -5.76756671071052551e-02  /**< LIO_4K x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_Y            -2.39010959863662720e-01   /**< LIO_4K y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_Z             9.69302475452423096e-01  /**< LIO_4K z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name Caltech 40m Prototype Detector constants
 * The following constants describe the location and geometry of the
 * Caltech 40m Prototype Detector.
 */
/*@{*/
#define LAL_CIT_40_DETECTOR_NAME                   "CIT_40"    /**< CIT_40 detector name string */
#define LAL_CIT_40_DETECTOR_PREFIX                 "P1"    /**< CIT_40 detector prefix string */
#define LAL_CIT_40_DETECTOR_LONGITUDE_RAD          -2.06175744538    /**< CIT_40 vertex longitude (rad) */
#define LAL_CIT_40_DETECTOR_LATITUDE_RAD           0.59637900541    /**< CIT_40 vertex latitude (rad) */
#define LAL_CIT_40_DETECTOR_ELEVATION_SI           0    /**< CIT_40 vertex elevation (m) */
#define LAL_CIT_40_DETECTOR_ARM_X_AZIMUTH_RAD      3.14159265359    /**< CIT_40 x arm azimuth (rad) */
#define LAL_CIT_40_DETECTOR_ARM_Y_AZIMUTH_RAD      1.57079632679    /**< CIT_40 y arm azimuth (rad) */
#define LAL_CIT_40_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< CIT_40 x arm altitude (rad) */
#define LAL_CIT_40_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< CIT_40 y arm altitude (rad) */
#define LAL_CIT_40_DETECTOR_ARM_X_MIDPOINT_SI      19.12500000000    /**< CIT_40 x arm midpoint (m) */
#define LAL_CIT_40_DETECTOR_ARM_Y_MIDPOINT_SI      19.12500000000    /**< CIT_40 y arm midpoint (m) */
#define LAL_CIT_40_VERTEX_LOCATION_X_SI            -2.49064958347e+06    /**< CIT_40 x-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_VERTEX_LOCATION_Y_SI            -4.65869968211e+06    /**< CIT_40 y-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_VERTEX_LOCATION_Z_SI            3.56206411403e+06    /**< CIT_40 z-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_ARM_X_DIRECTION_X               -0.26480331633    /**< CIT_40 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_X_DIRECTION_Y               -0.49530818538    /**< CIT_40 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_X_DIRECTION_Z               -0.82737476706    /**< CIT_40 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_X               0.88188012386    /**< CIT_40 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_Y               -0.47147369718    /**< CIT_40 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_Z               0.00000000000    /**< CIT_40 z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/**
 * \name Einstein Telescop 10km Interferometric Detector constants
 * The following constants describe the locations and geometrys of the
 * three 10km Interferometric Detectors for the planned third generation
 * Einstein Telescop detector as well as the theoretical null stream.
 * See T1400308
 */
/*@{*/
#define LAL_ET1_DETECTOR_NAME                      "ET1_T1400308"    /**< ET1 detector name string */
#define LAL_ET1_DETECTOR_PREFIX                    "E1"    /**< ET1 detector prefix string */
#define LAL_ET1_DETECTOR_LONGITUDE_RAD             0.18333805213    /**< ET1 vertex longitude (rad) */
#define LAL_ET1_DETECTOR_LATITUDE_RAD              0.76151183984    /**< ET1 vertex latitude (rad) */
#define LAL_ET1_DETECTOR_ELEVATION_SI              51.884    /**< ET1 vertex elevation (m) */
#define LAL_ET1_DETECTOR_ARM_X_AZIMUTH_RAD         0.33916285222    /**< ET1 x arm azimuth (rad) */
#define LAL_ET1_DETECTOR_ARM_Y_AZIMUTH_RAD         5.57515060820    /**< ET1 y arm azimuth (rad) */
#define LAL_ET1_DETECTOR_ARM_X_ALTITUDE_RAD        0.00000000000    /**< ET1 x arm altitude (rad) */
#define LAL_ET1_DETECTOR_ARM_Y_ALTITUDE_RAD        0.00000000000    /**< ET1 y arm altitude (rad) */
#define LAL_ET1_DETECTOR_ARM_X_MIDPOINT_SI         5000.00000000000    /**< ET1 x arm midpoint (m) */
#define LAL_ET1_DETECTOR_ARM_Y_MIDPOINT_SI         5000.00000000000    /**< ET1 y arm midpoint (m) */
#define LAL_ET1_VERTEX_LOCATION_X_SI               4.54637409900e+06    /**< ET1 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_VERTEX_LOCATION_Y_SI               8.42989697626e+05    /**< ET1 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_VERTEX_LOCATION_Z_SI               4.37857696241e+06    /**< ET1 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_ARM_X_DIRECTION_X                  -0.70045821479    /**< ET1 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_X_DIRECTION_Y                  0.20848948619    /**< ET1 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_X_DIRECTION_Z                  0.68256166277    /**< ET1 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_X                  -0.39681482542    /**< ET1 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_Y                  -0.73500471881    /**< ET1 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_Z                  0.54982366052    /**< ET1 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET2_DETECTOR_NAME                      "ET2_T1400308"    /**< ET2 detector name string */
#define LAL_ET2_DETECTOR_PREFIX                    "E2"    /**< ET2 detector prefix string */
#define LAL_ET2_DETECTOR_LONGITUDE_RAD             0.18405858870    /**< ET2 vertex longitude (rad) */
#define LAL_ET2_DETECTOR_LATITUDE_RAD              0.76299307990    /**< ET2 vertex latitude (rad) */
#define LAL_ET2_DETECTOR_ELEVATION_SI              59.735    /**< ET2 vertex elevation (m) */
#define LAL_ET2_DETECTOR_ARM_X_AZIMUTH_RAD         4.52795305701    /**< ET2 x arm azimuth (rad) */
#define LAL_ET2_DETECTOR_ARM_Y_AZIMUTH_RAD         3.48075550581    /**< ET2 y arm azimuth (rad) */
#define LAL_ET2_DETECTOR_ARM_X_ALTITUDE_RAD        0.00000000000    /**< ET2 x arm altitude (rad) */
#define LAL_ET2_DETECTOR_ARM_Y_ALTITUDE_RAD        0.00000000000    /**< ET2 y arm altitude (rad) */
#define LAL_ET2_DETECTOR_ARM_X_MIDPOINT_SI         5000.00000000000    /**< ET2 x arm midpoint (m) */
#define LAL_ET2_DETECTOR_ARM_Y_MIDPOINT_SI         5000.00000000000    /**< ET2 y arm midpoint (m) */
#define LAL_ET2_VERTEX_LOCATION_X_SI               4.53936951685e+06    /**< ET2 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_VERTEX_LOCATION_Y_SI               8.45074592488e+05    /**< ET2 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_VERTEX_LOCATION_Z_SI               4.38540257904e+06    /**< ET2 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_ARM_X_DIRECTION_X                  0.30364338937    /**< ET2 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_X_DIRECTION_Y                  -0.94349420500    /**< ET2 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_X_DIRECTION_Z                  -0.13273800225    /**< ET2 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_X                  0.70045821479    /**< ET2 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_Y                  -0.20848948619    /**< ET2 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_Z                  -0.68256166277    /**< ET2 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET3_DETECTOR_NAME                      "ET3_T1400308"    /**< ET3 detector name string */
#define LAL_ET3_DETECTOR_PREFIX                    "E3"    /**< ET3 detector prefix string */
#define LAL_ET3_DETECTOR_LONGITUDE_RAD             0.18192996730    /**< ET3 vertex longitude (rad) */
#define LAL_ET3_DETECTOR_LATITUDE_RAD              0.76270463257    /**< ET3 vertex latitude (rad) */
#define LAL_ET3_DETECTOR_ELEVATION_SI              59.727    /**< ET3 vertex elevation (m) */
#define LAL_ET3_DETECTOR_ARM_X_AZIMUTH_RAD         2.43355795462    /**< ET3 x arm azimuth (rad) */
#define LAL_ET3_DETECTOR_ARM_Y_AZIMUTH_RAD         1.38636040342    /**< ET3 y arm azimuth (rad) */
#define LAL_ET3_DETECTOR_ARM_X_ALTITUDE_RAD        0.00000000000    /**< ET3 x arm altitude (rad) */
#define LAL_ET3_DETECTOR_ARM_Y_ALTITUDE_RAD        0.00000000000    /**< ET3 y arm altitude (rad) */
#define LAL_ET3_DETECTOR_ARM_X_MIDPOINT_SI         5000.00000000000    /**< ET3 x arm midpoint (m) */
#define LAL_ET3_DETECTOR_ARM_Y_MIDPOINT_SI         5000.00000000000    /**< ET3 y arm midpoint (m) */
#define LAL_ET3_VERTEX_LOCATION_X_SI               4.54240595075e+06    /**< ET3 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_VERTEX_LOCATION_Y_SI               8.35639650438e+05    /**< ET3 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_VERTEX_LOCATION_Z_SI               4.38407519902e+06    /**< ET3 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_ARM_X_DIRECTION_X                  0.39681482542    /**< ET3 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_X_DIRECTION_Y                  0.73500471881    /**< ET3 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_X_DIRECTION_Z                  -0.54982366052    /**< ET3 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_X                  -0.30364338937    /**< ET3 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_Y                  0.94349420500    /**< ET3 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_Z                  0.13273800225    /**< ET3 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET0_DETECTOR_NAME                      "ET0_T1400308"    /**< ET0 detector name string */
#define LAL_ET0_DETECTOR_PREFIX                    "E0"    /**< ET0 detector prefix string */
#define LAL_ET0_DETECTOR_LONGITUDE_RAD             0.18192996730    /**< ET0 vertex longitude (rad) */
#define LAL_ET0_DETECTOR_LATITUDE_RAD              0.76270463257    /**< ET0 vertex latitude (rad) */
#define LAL_ET0_DETECTOR_ELEVATION_SI              59.727    /**< ET0 vertex elevation (m) */
#define LAL_ET0_DETECTOR_ARM_X_AZIMUTH_RAD         0.00000000000    /**< ET0 x arm azimuth (rad) */
#define LAL_ET0_DETECTOR_ARM_Y_AZIMUTH_RAD         0.00000000000    /**< ET0 y arm azimuth (rad) */
#define LAL_ET0_DETECTOR_ARM_X_ALTITUDE_RAD        0.00000000000    /**< ET0 x arm altitude (rad) */
#define LAL_ET0_DETECTOR_ARM_Y_ALTITUDE_RAD        0.00000000000    /**< ET0 y arm altitude (rad) */
#define LAL_ET0_DETECTOR_ARM_X_MIDPOINT_SI         0.00000000000    /**< ET0 x arm midpoint (m) */
#define LAL_ET0_DETECTOR_ARM_Y_MIDPOINT_SI         0.00000000000    /**< ET0 y arm midpoint (m) */
#define LAL_ET0_VERTEX_LOCATION_X_SI               4.54240595075e+06    /**< ET0 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_VERTEX_LOCATION_Y_SI               8.35639650438e+05    /**< ET0 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_VERTEX_LOCATION_Z_SI               4.38407519902e+06    /**< ET0 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_ARM_X_DIRECTION_X                  0.00000000000    /**< ET0 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_X_DIRECTION_Y                  0.00000000000    /**< ET0 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_X_DIRECTION_Z                  0.00000000000    /**< ET0 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_X                  0.00000000000    /**< ET0 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_Y                  0.00000000000    /**< ET0 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_Z                  0.00000000000    /**< ET0 z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/

/**
 * \name KAGRA Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * KAGRA Interferometric Detector.
 * \sa
 * > Yoshio Saito, "KAGRA location", KAGRA Technical Document JGW-G1503824
 * > http://gwdoc.icrr.u-tokyo.ac.jp/cgi-bin/DocDB/ShowDocument?docid=3824
 */
/*@{*/
#define LAL_KAGRA_DETECTOR_NAME                   "KAGRA"    /**< KAGRA detector name string */
#define LAL_KAGRA_DETECTOR_PREFIX                 "K1"    /**< KAGRA detector prefix string */
#define LAL_KAGRA_DETECTOR_LONGITUDE_RAD          2.396441015    /**< KAGRA vertex longitude (rad) */
#define LAL_KAGRA_DETECTOR_LATITUDE_RAD           0.6355068497    /**< KAGRA vertex latitude (rad) */
#define LAL_KAGRA_DETECTOR_ELEVATION_SI           414.181    /**< KAGRA vertex elevation (m) */
#define LAL_KAGRA_DETECTOR_ARM_X_AZIMUTH_RAD      1.054113    /**< KAGRA x arm azimuth (rad) */
#define LAL_KAGRA_DETECTOR_ARM_Y_AZIMUTH_RAD      -0.5166798    /**< KAGRA y arm azimuth (rad) */
#define LAL_KAGRA_DETECTOR_ARM_X_ALTITUDE_RAD     0.0031414    /**< KAGRA x arm altitude (rad) */
#define LAL_KAGRA_DETECTOR_ARM_Y_ALTITUDE_RAD     -0.0036270    /**< KAGRA y arm altitude (rad) */
#define LAL_KAGRA_DETECTOR_ARM_X_MIDPOINT_SI      1513.2535    /**< KAGRA x arm midpoint (m) */
#define LAL_KAGRA_DETECTOR_ARM_Y_MIDPOINT_SI      1511.611    /**< KAGRA y arm midpoint (m) */
#define LAL_KAGRA_VERTEX_LOCATION_X_SI            -3777336.024    /**< KAGRA x-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_VERTEX_LOCATION_Y_SI            3484898.411    /**< KAGRA y-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_VERTEX_LOCATION_Z_SI            3765313.697    /**< KAGRA z-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_ARM_X_DIRECTION_X               -0.3759040    /**< KAGRA x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_X_DIRECTION_Y               -0.8361583    /**< KAGRA y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_X_DIRECTION_Z               0.3994189    /**< KAGRA z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_X               0.7164378    /**< KAGRA x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_Y               0.01114076    /**< KAGRA y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_Z               0.6975620    /**< KAGRA z-component of unit vector pointing along y arm in Earth-centered frame */
/*@}*/


/* Resonant Mass (Bar) Detectors */


/**
 * \name ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis" constants
 * The following constants describe the location and geometry of the
 * ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis".
 */
/*@{*/
#define LAL_ALLEGRO_320_DETECTOR_NAME                   "ALLEGRO_320"    /**< ALLEGRO_320 detector name string */
#define LAL_ALLEGRO_320_DETECTOR_PREFIX                 "A1"    /**< ALLEGRO_320 detector prefix string */
#define LAL_ALLEGRO_320_DETECTOR_LONGITUDE_RAD          -1.59137068496    /**< ALLEGRO_320 vertex longitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_LATITUDE_RAD           0.53079879206    /**< ALLEGRO_320 vertex latitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ELEVATION_SI           0    /**< ALLEGRO_320 vertex elevation (m) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_AZIMUTH_RAD      -0.69813170080    /**< ALLEGRO_320 x arm azimuth (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_AZIMUTH_RAD      0.00000000000    /**< ALLEGRO_320 y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< ALLEGRO_320 x arm altitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< ALLEGRO_320 y arm altitude (rad) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_MIDPOINT_SI      0.00000000000    /**< ALLEGRO_320 x arm midpoint (m) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_MIDPOINT_SI      0.00000000000    /**< ALLEGRO_320 y arm midpoint (m) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_X_SI            -1.13258964140e+05    /**< ALLEGRO_320 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_Y_SI            -5.50408337391e+06    /**< ALLEGRO_320 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_Z_SI            3.20989567981e+06    /**< ALLEGRO_320 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_X                -0.63467362345    /**< ALLEGRO_320 x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_Y                0.40093077976    /**< ALLEGRO_320 y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_Z                0.66063901000    /**< ALLEGRO_320 z-component of unit vector pointing along axis in Earth-centered frame */
/*@}*/

/**
 * \name AURIGA Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * AURIGA Resonant Mass Detector.
 */
/*@{*/
#define LAL_AURIGA_DETECTOR_NAME                   "AURIGA"    /**< AURIGA detector name string */
#define LAL_AURIGA_DETECTOR_PREFIX                 "O1"    /**< AURIGA detector prefix string */
#define LAL_AURIGA_DETECTOR_LONGITUDE_RAD          0.20853775679    /**< AURIGA vertex longitude (rad) */
#define LAL_AURIGA_DETECTOR_LATITUDE_RAD           0.79156499342    /**< AURIGA vertex latitude (rad) */
#define LAL_AURIGA_DETECTOR_ELEVATION_SI           0    /**< AURIGA vertex elevation (m) */
#define LAL_AURIGA_DETECTOR_ARM_X_AZIMUTH_RAD      0.76794487088    /**< AURIGA x arm azimuth (rad) */
#define LAL_AURIGA_DETECTOR_ARM_Y_AZIMUTH_RAD      0.00000000000    /**< AURIGA y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< AURIGA x arm altitude (rad) */
#define LAL_AURIGA_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< AURIGA y arm altitude (rad) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_X_MIDPOINT_SI      0.00000000000    /**< AURIGA x arm midpoint (m) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_Y_MIDPOINT_SI      0.00000000000    /**< AURIGA y arm midpoint (m) UNUSED FOR BARS */
#define LAL_AURIGA_VERTEX_LOCATION_X_SI            4.39246733007e+06    /**< AURIGA x-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_VERTEX_LOCATION_Y_SI            9.29508666967e+05    /**< AURIGA y-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_VERTEX_LOCATION_Z_SI            4.51502913071e+06    /**< AURIGA z-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_AXIS_DIRECTION_X                -0.64450412225    /**< AURIGA x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_AURIGA_AXIS_DIRECTION_Y                0.57365538956    /**< AURIGA y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_AURIGA_AXIS_DIRECTION_Z                0.50550364038    /**< AURIGA z-component of unit vector pointing along axis in Earth-centered frame */
/*@}*/

/**
 * \name EXPLORER Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * EXPLORER Resonant Mass Detector.
 */
/*@{*/
#define LAL_EXPLORER_DETECTOR_NAME                   "EXPLORER"    /**< EXPLORER detector name string */
#define LAL_EXPLORER_DETECTOR_PREFIX                 "C1"            /**< EXPLORER detector prefix string */
#define LAL_EXPLORER_DETECTOR_LONGITUDE_RAD          0.10821041362    /**< EXPLORER vertex longitude (rad) */
#define LAL_EXPLORER_DETECTOR_LATITUDE_RAD           0.81070543755    /**< EXPLORER vertex latitude (rad) */
#define LAL_EXPLORER_DETECTOR_ELEVATION_SI           0    /**< EXPLORER vertex elevation (m) */
#define LAL_EXPLORER_DETECTOR_ARM_X_AZIMUTH_RAD      0.68067840828    /**< EXPLORER x arm azimuth (rad) */
#define LAL_EXPLORER_DETECTOR_ARM_Y_AZIMUTH_RAD      0.00000000000    /**< EXPLORER y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< EXPLORER x arm altitude (rad) */
#define LAL_EXPLORER_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< EXPLORER y arm altitude (rad) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_X_MIDPOINT_SI      0.00000000000    /**< EXPLORER x arm midpoint (m) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_Y_MIDPOINT_SI      0.00000000000    /**< EXPLORER y arm midpoint (m) UNUSED FOR BARS */
#define LAL_EXPLORER_VERTEX_LOCATION_X_SI            4.37645395452e+06    /**< EXPLORER x-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_VERTEX_LOCATION_Y_SI            4.75435044067e+05    /**< EXPLORER y-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_VERTEX_LOCATION_Z_SI            4.59985274450e+06    /**< EXPLORER z-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_AXIS_DIRECTION_X                -0.62792641437    /**< EXPLORER x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_EXPLORER_AXIS_DIRECTION_Y                0.56480832712    /**< EXPLORER y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_EXPLORER_AXIS_DIRECTION_Z                0.53544371484    /**< EXPLORER z-component of unit vector pointing along axis in Earth-centered frame */
/*@}*/

/**
 * \name Nautilus Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * Nautilus Resonant Mass Detector.
 */
/*@{*/
#define LAL_NAUTILUS_DETECTOR_NAME                   "Nautilus"    /**< Nautilus detector name string */
#define LAL_NAUTILUS_DETECTOR_PREFIX                 "N1"    /**< Nautilus detector prefix string */
#define LAL_NAUTILUS_DETECTOR_LONGITUDE_RAD          0.22117684946    /**< Nautilus vertex longitude (rad) */
#define LAL_NAUTILUS_DETECTOR_LATITUDE_RAD           0.72996456710    /**< Nautilus vertex latitude (rad) */
#define LAL_NAUTILUS_DETECTOR_ELEVATION_SI           0    /**< Nautilus vertex elevation (m) */
#define LAL_NAUTILUS_DETECTOR_ARM_X_AZIMUTH_RAD      0.76794487088    /**< Nautilus x arm azimuth (rad) */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_AZIMUTH_RAD      0.00000000000    /**< Nautilus y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< Nautilus x arm altitude (rad) */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< Nautilus y arm altitude (rad) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_X_MIDPOINT_SI      0.00000000000    /**< Nautilus x arm midpoint (m) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_MIDPOINT_SI      0.00000000000    /**< Nautilus y arm midpoint (m) UNUSED FOR BARS */
#define LAL_NAUTILUS_VERTEX_LOCATION_X_SI            4.64410999868e+06    /**< Nautilus x-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_VERTEX_LOCATION_Y_SI            1.04425342477e+06    /**< Nautilus y-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_VERTEX_LOCATION_Z_SI            4.23104713307e+06    /**< Nautilus z-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_AXIS_DIRECTION_X                -0.62039441384    /**< Nautilus x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NAUTILUS_AXIS_DIRECTION_Y                0.57250373141    /**< Nautilus y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NAUTILUS_AXIS_DIRECTION_Z                0.53605060283    /**< Nautilus z-component of unit vector pointing along axis in Earth-centered frame */
/*@}*/

/**
 * \name NIOBE Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * NIOBE Resonant Mass Detector.
 */
/*@{*/
#define LAL_NIOBE_DETECTOR_NAME                   "NIOBE"    /**< NIOBE detector name string */
#define LAL_NIOBE_DETECTOR_PREFIX                 "B1"    /**< NIOBE detector prefix string */
#define LAL_NIOBE_DETECTOR_LONGITUDE_RAD          2.02138216202    /**< NIOBE vertex longitude (rad) */
#define LAL_NIOBE_DETECTOR_LATITUDE_RAD           -0.55734180780    /**< NIOBE vertex latitude (rad) */
#define LAL_NIOBE_DETECTOR_ELEVATION_SI           0    /**< NIOBE vertex elevation (m) */
#define LAL_NIOBE_DETECTOR_ARM_X_AZIMUTH_RAD      0.00000000000    /**< NIOBE x arm azimuth (rad) */
#define LAL_NIOBE_DETECTOR_ARM_Y_AZIMUTH_RAD      0.00000000000    /**< NIOBE y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_X_ALTITUDE_RAD     0.00000000000    /**< NIOBE x arm altitude (rad) */
#define LAL_NIOBE_DETECTOR_ARM_Y_ALTITUDE_RAD     0.00000000000    /**< NIOBE y arm altitude (rad) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_X_MIDPOINT_SI      0.00000000000    /**< NIOBE x arm midpoint (m) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_Y_MIDPOINT_SI      0.00000000000    /**< NIOBE y arm midpoint (m) UNUSED FOR BARS */
#define LAL_NIOBE_VERTEX_LOCATION_X_SI            -2.35948871453e+06    /**< NIOBE x-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_VERTEX_LOCATION_Y_SI            4.87721571259e+06    /**< NIOBE y-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_VERTEX_LOCATION_Z_SI            -3.35416003274e+06    /**< NIOBE z-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_AXIS_DIRECTION_X                -0.23034623759    /**< NIOBE x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NIOBE_AXIS_DIRECTION_Y                0.47614056486    /**< NIOBE y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NIOBE_AXIS_DIRECTION_Z                0.84866411101    /**< NIOBE z-component of unit vector pointing along axis in Earth-centered frame */
/*@}*/


/* LALDetector Macro */
#define LAL_CAT(x,y) x ## y
#define LAL_XCAT(x,y) LAL_CAT(x,y)

/** expands to constant c of detector d */
#define LAL_DETECTOR_CONSTANT(d,c) LAL_XCAT(LAL_XCAT(LAL_,d),LAL_XCAT(_,c))

/** initializer for detector location vector */
#define LAL_DETECTOR_LOCATION(d) \
{ \
LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_X_SI),\
LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Y_SI),\
LAL_DETECTOR_CONSTANT(d,VERTEX_LOCATION_Z_SI) \
}

/** expands to component c (X,Y,Z) of arm X of detector d */
#define LAL_ARM_X(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_X_DIRECTION_,c))

/** expands to component c (X,Y,Z) of arm Y of detector d */
#define LAL_ARM_Y(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(ARM_Y_DIRECTION_,c))

/** expands to component c (X,Y,Z) of axis of detector d */
#define LAL_AXIS(d,c) LAL_DETECTOR_CONSTANT(d,LAL_XCAT(AXIS_DIRECTION_,c))

/** expands to a 3x3 matix initializer for the response for IFODIFF detector d */
#define LAL_DETECTOR_RESPONSE_IFODIFF(d) \
{ \
{ \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z) )  \
}, \
{ \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z) )  \
}, \
{ \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,X) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z) - LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z) )  \
} \
}

/** expands to a 3x3 matix initializer for the response for IFOCOMM detector d */
#define LAL_DETECTOR_RESPONSE_IFOCOMM(d) \
{ \
{ \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,X) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z) )  \
}, \
{ \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z) )  \
}, \
{ \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,X) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X) ), \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y) ), \
0.5*( LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z) + LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z) )  \
} \
}

/** expands to a 3x3 matix initializer for the response for IFOXARM detector d */
#define LAL_DETECTOR_RESPONSE_IFOXARM(d) \
{ \
{ \
0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,X), \
0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,Y), \
0.5 * LAL_ARM_X(d,X) * LAL_ARM_X(d,Z)  \
}, \
{ \
0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,X), \
0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,Y), \
0.5 * LAL_ARM_X(d,Y) * LAL_ARM_X(d,Z)  \
}, \
{ \
0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,X), \
0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,Y), \
0.5 * LAL_ARM_X(d,Z) * LAL_ARM_X(d,Z)  \
} \
}

/** expands to a 3x3 matix initializer for the response for IFOYARM detector d */
#define LAL_DETECTOR_RESPONSE_IFOYARM(d) \
{ \
{ \
0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,X), \
0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Y), \
0.5 * LAL_ARM_Y(d,X) * LAL_ARM_Y(d,Z)  \
}, \
{ \
0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,X), \
0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Y), \
0.5 * LAL_ARM_Y(d,Y) * LAL_ARM_Y(d,Z)  \
}, \
{ \
0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,X), \
0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Y), \
0.5 * LAL_ARM_Y(d,Z) * LAL_ARM_Y(d,Z)  \
} \
}

/** expands to a 3x3 matix initializer for the response for CYLBAR detector d */
#define LAL_DETECTOR_RESPONSE_CYLBAR(d) \
{ \
{ \
LAL_AXIS(d,X) * LAL_AXIS(d,X), \
LAL_AXIS(d,X) * LAL_AXIS(d,Y), \
LAL_AXIS(d,X) * LAL_AXIS(d,Z)  \
}, \
{ \
LAL_AXIS(d,Y) * LAL_AXIS(d,X), \
LAL_AXIS(d,Y) * LAL_AXIS(d,Y), \
LAL_AXIS(d,Y) * LAL_AXIS(d,Z)  \
}, \
{ \
LAL_AXIS(d,Z) * LAL_AXIS(d,X), \
LAL_AXIS(d,Z) * LAL_AXIS(d,Y), \
LAL_AXIS(d,Z) * LAL_AXIS(d,Z)  \
} \
}

#define LAL_FR_STREAM_DETECTOR_STRUCT(d) \
{ \
LAL_DETECTOR_CONSTANT(d,DETECTOR_NAME), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_PREFIX), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_LONGITUDE_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_LATITUDE_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ELEVATION_SI), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_ALTITUDE_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_AZIMUTH_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_ALTITUDE_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_AZIMUTH_RAD), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_X_MIDPOINT_SI), \
LAL_DETECTOR_CONSTANT(d,DETECTOR_ARM_Y_MIDPOINT_SI) \
}

#define LAL_DETECTOR_RESPONSE(d,t) \
LAL_XCAT( LAL_DETECTOR_RESPONSE_, t )(d)

#define LAL_DETECTOR_STRUCT(d,t) \
{ \
LAL_DETECTOR_LOCATION(d),      \
LAL_DETECTOR_RESPONSE(d,t),    \
LAL_XCAT(LALDETECTORTYPE_,t),  \
LAL_FR_STREAM_DETECTOR_STRUCT(d)      \
}

typedef enum tagLALDetectorType {
    LALDETECTORTYPE_ABSENT,        /**< No FrDetector associated with this detector */
    LALDETECTORTYPE_IFODIFF,    /**< IFO in differential mode */
    LALDETECTORTYPE_IFOXARM,    /**< IFO in one-armed mode (X arm) */
    LALDETECTORTYPE_IFOYARM,    /**< IFO in one-armed mode (Y arm) */
    LALDETECTORTYPE_IFOCOMM,    /**< IFO in common mode */
    LALDETECTORTYPE_CYLBAR        /**< Cylindrical bar */
}
LALDetectorType;




typedef struct tagLALFrDetector
{
    char    name[LALNameLength];    /**< A unique identifying string */
    char    prefix[3];        /**< Two-letter prefix for detector's channel names */
    double    vertexLongitudeRadians;    /**< The geodetic longitude \f$\lambda\f$ of the vertex in radians */
    double    vertexLatitudeRadians;    /**< The geodetic latitude \f$\beta\f$ of the vertex in radians */
    float    vertexElevation;    /**< The height of the vertex above the reference ellipsoid in meters */
    float    xArmAltitudeRadians;    /**< The angle \f${\mathcal{A}}_X\f$ up from the local tangent plane of the reference ellipsoid to the X arm (or bar's cylidrical axis) in radians */
    float    xArmAzimuthRadians;    /**< The angle \f$\zeta_X\f$ clockwise from North to the projection of the X arm (or bar's cylidrical axis) into the local tangent plane of the reference ellipsoid in radians */
    float    yArmAltitudeRadians;    /**< The angle \f${\mathcal{A}}_Y\f$ up from the local tangent plane of the reference ellipsoid to the Y arm in radians (unused for bars: set it to zero) */
    float    yArmAzimuthRadians;    /**< The angle \f$\zeta_Y\f$ clockwise from North to the projection of the Y arm into the local tangent plane of the reference ellipsoid in radians (unused for bars: set it to zero) */
    float    xArmMidpoint;        /**< The distance to the midpoint of the X arm in meters (unused for bars: set it to zero) */
    float    yArmMidpoint;        /**< The distance to the midpoint of the Y arm in meters (unused for bars: set it to zero) */
}
LALFrDetector;

typedef struct tagLALDetector
{
    double        location[3];    /**< The three components, in an Earth-fixed Cartesian coordinate system, of the position vector from the center of the Earth to the detector in meters */
    double        response[3][3];    /**< The Earth-fixed Cartesian components of the detector's response tensor \f$d^{ab}\f$ */
    LALDetectorType    type;        /**< The type of the detector (e.g., IFO in differential mode, cylindrical bar, etc.) */
    LALFrDetector    frDetector;    /**< The original LALFrDetector structure from which this was created */
}
LALDetector;

/* Functions */
void antenna_pattern_comm(const double    D[3][3],
                          const double    ra,
                          const double    dec,
                          const double    psi,
                          const double    gmst,
                          double    *fplus,
                          double    *fcross);

double time_delay_from_earth_center(const double    detloc[3],
                                    const double          ra,
                                    const double          de,
                                    const double    gpstime);

void set_detector(const char *prefix, const LALDetector **detector);

#endif

