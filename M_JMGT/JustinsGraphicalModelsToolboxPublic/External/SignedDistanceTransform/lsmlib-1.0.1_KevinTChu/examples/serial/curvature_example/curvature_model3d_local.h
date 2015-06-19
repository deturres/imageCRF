/*
 * File:        curvature_model3d_local.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 151 $
 * Modified:    $Date: 2009-01-27 09:03:59 -0800 (Tue, 27 Jan 2009) $
 * Description: Support header file for localized 3D constant curvature flow.
 */
#ifndef INCLUDED_CURV_MODEL3D_LOCAL_H
#define INCLUDED_CURV_MODEL3D_LOCAL_H

void  curvatureModelMedium3dLocalMainLoop(Options *,LSM_DataArrays  *,Grid  *,FILE *);
void  reinitializeMedium3dLocal(LSM_DataArrays *,Grid *,Options *,LSMLIB_REAL);

#endif
