/*
 * File:        lsm_FMM_field_extension2d.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Implementation of 2D Fast Marching Method for computing
 *              signed distance functions and extension fields
 */
 

/*
 * lsm_FMM_field_extension2d.c makes use of the generic implementation 
 * of the Fast Marching Method algorithm for computing signed distance
 * functions and extension fields provided by lsm_FMM_field_extension.c.
 */

#include "lsm_fast_marching_method.h"


/* Define required macros */
#define FMM_NDIM                         2
#define FMM_COMPUTE_DISTANCE_FUNCTION    computeDistanceFunction2d
#define FMM_COMPUTE_EXTENSION_FIELDS     computeExtensionFields2d
#define FMM_INITIALIZE_FRONT_ORDER1                                         \
        FMM_initializeFront_FieldExtension2d_Order1
#define FMM_INITIALIZE_FRONT_ORDER2                                         \
        FMM_initializeFront_FieldExtension2d_Order2
#define FMM_UPDATE_GRID_POINT_ORDER1                                        \
        FMM_updateGridPoint_FieldExtension2d_Order1
#define FMM_UPDATE_GRID_POINT_ORDER2                                        \
        FMM_updateGridPoint_FieldExtension2d_Order2


/* Include "templated" implementation of Fast Marching Method */
/* signed distance functions and extension fields.            */
#include "lsm_FMM_field_extension.c"

