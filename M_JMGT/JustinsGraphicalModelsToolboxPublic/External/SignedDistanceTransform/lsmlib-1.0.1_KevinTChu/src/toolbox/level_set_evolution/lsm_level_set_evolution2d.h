/*
 * File:        lsm_level_set_evolution2d.h
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Header file for 2D Fortran 77 level set evolution equation
 *              subroutines
 */

#ifndef INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_H
#define INCLUDED_LSM_LEVEL_SET_EVOLUTION_2D_H

#include "LSMLIB_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file lsm_level_set_evolution2d.h
 *
 * \brief
 * @ref lsm_level_set_evolution2d.h provides support for contributing to the
 * right-hand side of the level set evolution equation in two space
 * dimensions.
 * 
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                                name in
 *      C/C++ code                             Fortran code
 *      ----------                             ------------
 */
#define LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS       lsm2dzerooutlevelseteqnrhs_
#define LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS    lsm2daddadvectiontermtolserhs_
#define LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS   lsm2daddnormalveltermtolserhs_
#define LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS   \
                                          lsm2daddconstnormalveltermtolserhs_
#define LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS   lsm2daddconstcurvtermtolserhs_
#define LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS \
                                     lsm2daddconstprecomputedcurvtermtolserhs_	
#define LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS \
                                    lsm2daddexternalandnormalveltermtolserhs_				

/*!
 * LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS() zeros out the right-hand side of 
 * the level set equation when it is written in the form:
 *
 * \f[
 *
 *  \phi_t = ...
 *
 * \f]
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - *_gb (in):         index range for ghostbox
 *
 * Return value:         none
 */
void LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb);


/*!
 * LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS() adds the contribution of an 
 * advection term (external vector velocity field) to the right-hand 
 * side of the level set equation when it is written in the form:
 *
 * \f[
 *
 *    \phi_t = -\vec{V} \cdot \nabla \phi + ...
 *
 * \f]
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_* (in):        components of \f$ \nabla \phi\f$ at t = t_cur
 *  - vel_* (in):        components of velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb, 
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb, 
  const int *jhi_grad_phi_gb,
  const LSMLIB_REAL *vel_x, 
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb);


/*!
 * LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution of a 
 * normal (scalar) velocity term to the right-hand side of the level 
 * set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -V_n |\nabla \phi| + ... 
 *   
 * \f]
 *   
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_*_plus (in):   components of forward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - phi_*_minus (in):  components of backward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - vel_n (in):        normal velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb);


/*!
 * LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution 
 * of a normal (scalar) velocity term to the right-hand side of the 
 * level set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -V_n |\nabla \phi| + ... 
 *   
 * \f]
 *  
 * where the normal velocity, \f$ V_n \f$ is a constant.
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_*_plus (in):   components of forward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - phi_*_minus (in):  components of backward approx to \f$ \nabla \phi \f$
 *                       at t = t_cur
 *  - vel_n (in):        constant normal velocity at t = t_cur
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb);
   

/*!
 * LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS() adds the contribution 
 * of a normal (scalar) velocity term to the right-hand side of the 
 * level set equation when it is written in the form:
 *   
 * \f[
 *   
 *    \phi_t = -b kappa |\nabla \phi| + ... 
 *   
 * \f]
 *  
 * where the \f$ kappa \f$ is the mean curvature and \f$ b \f$ is a 
 * constant.
 *
 * Arguments:
 *  - lse_rhs (in/out):  right-hand of level set equation
 *  - phi_* (in):        first- and second-order partial derivatives 
 *                       of \f$ \phi \f$ 
 *  - b (in):            proporationality constant relating curvature
 *                       to the normal velocity
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void   LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS(
  const LSMLIB_REAL  *lse_rhs,
  const int *ilo_lse_rhs_gb, 
  const int *ihi_lse_rhs_gb,
  const int *jlo_lse_rhs_gb, 
  const int *jhi_lse_rhs_gb,
  const LSMLIB_REAL *phi_x, 
  const LSMLIB_REAL *phi_y,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const  LSMLIB_REAL *phi_xx,
  const  LSMLIB_REAL *phi_xy,
  const  LSMLIB_REAL *phi_yy,
  const  int *ilo_grad2_phi_gb,
  const  int *ihi_grad2_phi_gb,
  const  int *jlo_grad2_phi_gb, 
  const  int *jhi_grad2_phi_gb,
  const LSMLIB_REAL *b,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb);  

/*!
*
*  LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS() adds the contribution of 
*  a curvature term to the right-hand side of the level set equation 
*  when it is  written in the form:
*
* \f[
*   
*    \phi_t = -b kappa |\nabla \phi| + ... 
*   
* \f]
*  Arguments:
*    lse_rhs (in/out):  right-hand of level set equation
*    kappa (in):       curvature data array (assumed precomputed)
*    grad_phi_mag(in): gradient magnitude data_array
*    b     (in):        scalar curvature term component 
*    *_gb (in):         index range for ghostbox
*    *_fb (in):         index range for fillbox
*/

 void LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS(
  const LSMLIB_REAL *lse_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  LSMLIB_REAL *kappa,
  const int *ilo_kappa_gb, 
  const int *ihi_kappa_gb,
  const int *jlo_kappa_gb, 
  const int *jhi_kappa_gb,
  const LSMLIB_REAL *grad_mag_phi,
  const int *ilo_phi_gb, 
  const int *ihi_phi_gb,
  const int *jlo_phi_gb, 
  const int *jhi_phi_gb, 
  const LSMLIB_REAL *b, 
  const int *ilo_rhs_fb,
  const int *ihi_rhs_fb,
  const int *jlo_rhs_fb,
  const int *jhi_rhs_fb);
   
/*!
*
*  LSERHSLSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS() adds the contribution 
*  of a normal (scalar) velocity term AND advective term to the right-hand 
*  side of the level set equation when it is written in the form:
*
*   \f[
*    \phi_t = (- V  -  V_n \nabla \phi/|\nabla \phi|) . \nabla \phi
*   i.e.
*   \phi_t = - V . \nabla \phi - V_n |\nabla \phi|
*
*  \f]
*
*  Note that the upwinding selection is made assuming \f$ \nabla \phi \f$ is 1 
*  which should be a reasonable choice if \f$ \phi \f$ is close to a signed 
*  distance function.
*  See discussion 'Adding an External Velocity field' on Pg 59, Osher/Fedkiw 
*  book.
*
*  Arguments:
*    lse_rhs (in/out):  right-hand of level set equation
*    phi_*_plus (in):   components of forward approx to \f$ \nabla \phi \f$ at 
*                       t = t_cur
*    phi_*_minus (in):  components of backward approx to \f$ \nabla \phi \f$at 
*                       t = t_cur
*    vel_n (in):        normal velocity at t = t_cur
*    vel_[xy]:          external velocity (advective), V = (vel_x, vel_y)
*    *_gb (in):         index range for ghostbox
*    *_fb (in):         index range for fillbox
*
*/  
void LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS(
  LSMLIB_REAL *lse_rhs,
  const int *ilo_rhs_gb, 
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb, 
  const int *jhi_rhs_gb,
  const LSMLIB_REAL *phi_x_plus, 
  const LSMLIB_REAL *phi_y_plus,
  const int *ilo_grad_phi_plus_gb, 
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb, 
  const int *jhi_grad_phi_plus_gb,
  const LSMLIB_REAL *phi_x_minus, 
  const LSMLIB_REAL *phi_y_minus,
  const int *ilo_grad_phi_minus_gb, 
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb, 
  const int *jhi_grad_phi_minus_gb,
  const LSMLIB_REAL *vel_n,
  const LSMLIB_REAL *vel_x,
  const LSMLIB_REAL *vel_y,
  const int *ilo_vel_gb, 
  const int *ihi_vel_gb,
  const int *jlo_vel_gb, 
  const int *jhi_vel_gb,
  const int *ilo_rhs_fb, 
  const int *ihi_rhs_fb,
  const int *jlo_rhs_fb, 
  const int *jhi_rhs_fb);
#ifdef __cplusplus
}
#endif

#endif
