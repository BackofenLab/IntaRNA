#ifndef VIENNA_RNA_PACKAGE_CONFIG_H
#define VIENNA_RNA_PACKAGE_CONFIG_H

/* version number */
#define VRNA_VERSION  "2.6.4"

#define VRNA_VERSION_MAJOR  2
#define VRNA_VERSION_MINOR  6
#define VRNA_VERSION_PATCH  4

/*
 * The following pre-processor definitions specify whether
 * or not certain features were activated upon build-time
 */

/*
 * Build with deactivated C11 Features
 *
 * It this feature is missing, the next line defines
 * 'VRNA_DISABLE_C11_FEATURES'
 */


/*
 * Build with OpenMP support
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_OPENMP'
 */
#define VRNA_WITH_OPENMP

/*
 * Build with single precision partition function
 *
 * If this feature is present, the next line defines
 * 'USE_FLOAT_PF'
 */


/*
 * Build with JSON input/output support
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_JSON_SUPPORT'
 *
 * Note: ViennaRNA >= 2.6.3 is always built with JSON support.
 * This feature is still present to keep third-party code relying on it working.
 */
#define VRNA_WITH_JSON_SUPPORT

/*
 * Build with Support Vector Machine (SVM) Z-score feature in RNALfold
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_SVM'
 */
#define VRNA_WITH_SVM

/*
 * Build with GSL minimizers
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_GSL'
 */
#define VRNA_WITH_GSL

/*
 * Build with colored TTY output
 *
 * If this feature is missing, the next line defines
 * 'VRNA_WITHOUT_TTY_COLORS'
 */


/*
 * Build with Link Time Optimization support
 *
 * If this feature is enabled, the next line defines
 * 'VRNA_WITH_LTO'
 */


/*
 * Build with Naview Layout algorithm of Bruccoleri 1988
 *
 * If this feature is enabled, the next line defines
 * 'VRNA_WITH_NAVIEW_LAYOUT'
 */
#define VRNA_WITH_NAVIEW_LAYOUT

#endif
