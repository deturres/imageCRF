/*
 * File:        LevelSetFunctionIntegrator.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Implementation file for level set method integrator class
 */

#ifndef included_LevelSetFunctionIntegrator_cc
#define included_LevelSetFunctionIntegrator_cc

// System Headers
#include <float.h>
#include <sstream>

#include "LevelSetFunctionIntegrator.h" 
#include "LevelSetMethodPatchStrategy.h" 
#include "LevelSetMethodVelocityFieldStrategy.h" 
#include "LSMLIB_DefaultParameters.h" 
#include "BoundaryConditionModule.h" 

// SAMRAI header files
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "PatchLevel.h"
#include "RefineOperator.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

// headers for level set method numerical kernels
extern "C" {
  #include "lsm_level_set_evolution1d.h"
  #include "lsm_level_set_evolution2d.h"
  #include "lsm_level_set_evolution3d.h"
  #include "lsm_utilities1d.h"
  #include "lsm_utilities2d.h"
  #include "lsm_utilities3d.h"
}

// SAMRAI namespaces
using namespace pdat;

// VERSION INFORMATION
#define LEVEL_SET_METHOD_INTEGRATOR_VERSION	(1)


/*
 * Default values for user-defined parameters
 */
#define LSM_DEFAULT_START_TIME                           (0.0)
#define LSM_DEFAULT_END_TIME                             (10.0)
#define LSM_DEFAULT_REINITIALIZATION_INTERVAL            (10)
#define LSM_DEFAULT_REINITIALIZATION_MAX_ITERS           (25)
#define LSM_DEFAULT_ORTHOGONALIZATION_INTERVAL           (10)
#define LSM_DEFAULT_ORTHOGONALIZATION_MAX_ITERS          (25)
#define LSM_DEFAULT_USE_AMR                              (false)
#define LSM_DEFAULT_REGRID_INTERVAL                      (5)  // KTC - ADJUST
#define LSM_DEFAULT_TAG_BUFFER_WIDTH                     (2)  // KTC - ADJUST
#define LSM_DEFAULT_REFINEMENT_CUTOFF_VALUE              (1.0)  // KTC - ADJUST
#define LSM_DEFAULT_VERBOSE_MODE                         (false)
#define LSM_STOP_TOLERANCE_MAX_ITERATIONS                (1000)

#ifdef LSMLIB_DEBUG_NO_INLINE
#include "LevelSetFunctionIntegrator.inline"
#endif

/****************************************************************
 *
 * Implementation for LevelSetFunctionIntegrator methods.
 *
 ****************************************************************/

namespace LSMLIB {

/* Constructor */
template <int DIM> 
LevelSetFunctionIntegrator<DIM>::LevelSetFunctionIntegrator(
  Pointer<Database> input_db,
  Pointer< PatchHierarchy<DIM> > patch_hierarchy,
  LevelSetMethodPatchStrategy<DIM>* lsm_patch_strategy,
  LevelSetMethodVelocityFieldStrategy<DIM>* lsm_velocity_field_strategy,
  const int num_level_set_fcn_components,
  const int codimension,
  const string& object_name) 
:
  d_grad_psi_plus_handle(-1),
  d_grad_psi_minus_handle(-1),
  d_grad_psi_upwind_handle(-1),
  d_rhs_psi_handle(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!input_db.isNull());
  assert(!patch_hierarchy.isNull());
  assert(lsm_patch_strategy);
  assert(lsm_velocity_field_strategy);
  assert(!object_name.empty());
#endif

  // Check to make sure that DIM is 1, 2, or 3.
  // Dimensions greater than 3 are currently unsupported.
  if (DIM > 3) {
    TBOX_ERROR(  d_object_name 
              << "::LevelSetFunctionIntegrator(): "
              << "DIM > 3 not supported."
              << endl );
  } 
  if (DIM < 1) {
    TBOX_ERROR(  d_object_name 
              << "::LevelSetFunctionIntegrator(): "
              << "DIM must be positive."
              << endl );
  }

  // set pointers to major objects 
  d_patch_hierarchy = patch_hierarchy; 
  d_grid_geometry = d_patch_hierarchy->getGridGeometry();
  d_lsm_patch_strategy = lsm_patch_strategy;
  d_lsm_velocity_field_strategy = lsm_velocity_field_strategy;
  d_object_name = object_name;

  // set number of level set function components
  d_num_level_set_fcn_components = num_level_set_fcn_components;

  // set codimension for problem
  d_codimension = codimension;

  // create empty BoundaryConditionModule
  d_bc_module = new BoundaryConditionModule<DIM>;

  // initialize boundary condition data
  d_lower_bc_phi.resizeArray(d_num_level_set_fcn_components);
  d_upper_bc_phi.resizeArray(d_num_level_set_fcn_components);
  for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
    for (int dim=0; dim<DIM; dim++) {
      d_lower_bc_phi[comp](dim) = BoundaryConditionModule<DIM>::NONE;
      d_upper_bc_phi[comp](dim) = BoundaryConditionModule<DIM>::NONE;
    }
  }
  
  if (d_codimension == 2) {
    d_lower_bc_psi.resizeArray(d_num_level_set_fcn_components);
    d_upper_bc_psi.resizeArray(d_num_level_set_fcn_components);
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      for (int dim=0; dim<DIM; dim++) {
        d_lower_bc_psi[comp](dim) = BoundaryConditionModule<DIM>::NONE;
        d_upper_bc_psi[comp](dim) = BoundaryConditionModule<DIM>::NONE;
      }
    }
  }

  // register LevelSetFunctionIntegrator for restart
  RestartManager::getManager()->registerRestartItem(d_object_name, this);

  // initialize user-defined parameters from given input & restart databases
  bool is_from_restart = RestartManager::getManager()->isFromRestart();
  if (is_from_restart){
    getFromRestart();
  }
  getFromInput(input_db, is_from_restart);

  // initialize current time, integrator step and counter variables
  if (!is_from_restart) {
    d_current_time = d_start_time;
    d_num_integration_steps_taken = 0;
    d_regrid_count = 0;
    d_reinitialization_count = 0;
    d_orthogonalization_count = 0;
  }

  // initialize variables and communication objects
  initializeVariables();
  initializeCommunicationObjects();

  // create reinitialization algorithm for phi
  d_phi_reinitialization_alg = 
    new ReinitializationAlgorithm<DIM>(
      d_patch_hierarchy,
      d_phi_handles[0],
      d_control_volume_handle,
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      d_tvd_runge_kutta_order,
      d_cfl_number,
      d_reinitialization_stop_dist,
      d_reinitialization_max_iters,
      d_reinitialization_stop_tol,
      d_verbose_mode,
      "phi reinitialization algorithm");

  // create reinitialization algorithm for psi (if necessary)
  if (d_codimension == 2) {
    d_psi_reinitialization_alg = 
      new ReinitializationAlgorithm<DIM>(
        d_patch_hierarchy,
        d_psi_handles[0],
        d_control_volume_handle,
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        d_tvd_runge_kutta_order,
        d_cfl_number,
        d_reinitialization_stop_dist,
        d_reinitialization_max_iters,
        d_reinitialization_stop_tol,
        d_verbose_mode,
        "psi reinitialization algorithm");
   }
 
  // create orthogonalization algorithm for codimension-two problems
  if (d_codimension == 2) {
    d_orthogonalization_alg = 
      new OrthogonalizationAlgorithm<DIM>(
        d_patch_hierarchy,
        d_phi_handles[0],
        d_psi_handles[0],
        d_control_volume_handle,
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        d_tvd_runge_kutta_order,
        d_cfl_number,
        d_orthogonalization_stop_dist,
        d_orthogonalization_max_iters,
        d_orthogonalization_stop_tol,
        d_verbose_mode,
        "orthogonalization algorithm");
  } else {
    d_orthogonalization_alg.setNull(); 
  }
  d_orthogonalization_evolved_field = PHI;  // first orthogonalization
                                            // evolves phi

}


/* Destructor */
template <int DIM> 
LevelSetFunctionIntegrator<DIM>::~LevelSetFunctionIntegrator()
{
  // unregister this object as a restart item
  RestartManager::getManager()->unregisterRestartItem(d_object_name);

  // deallocate solution variables and persistent variables
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData(d_solution_variables); 
    level->deallocatePatchData(d_persistent_variables); 
  }

}


/* printClassData() */
template <int DIM> void 
LevelSetFunctionIntegrator<DIM>::printClassData(
  ostream& os) const
{
  os << "\n===================================" << endl;
  os << "LevelSetFunctionIntegrator<DIM>:" << endl;
  os << "(LevelSetFunctionIntegrator*) this = " 
     << (LevelSetFunctionIntegrator*)this << endl;
  os << "d_object_name = " << d_object_name << endl;

  os << "Level set method parameters" << endl;
  os << "---------------------------" << endl;
  os << "d_codimension = " << d_codimension << endl;
  os << "d_start_time = " << d_start_time << endl;
  os << "d_end_time = " << d_end_time << endl;
  os << "d_cfl_number = " << d_cfl_number << endl;
  os << "d_spatial_derivative_type = " << d_spatial_derivative_type << endl;
  os << "d_spatial_derivative_order = " << d_spatial_derivative_order << endl;
  os << "d_tvd_runge_kutta_order = " << d_tvd_runge_kutta_order << endl;
  os << "d_reinitialization_interval = " 
     << d_reinitialization_interval << endl;
  os << "d_reinitialization_stop_tol = " 
     << d_reinitialization_stop_tol << endl;
  os << "d_reinitialization_stop_dist = " 
     << d_reinitialization_stop_dist << endl;
  os << "d_reinitialization_max_iters = " 
     << d_reinitialization_max_iters << endl;
  os << "d_orthogonalization_interval = " 
     << d_orthogonalization_interval << endl;
  os << "d_orthogonalization_stop_tol = " 
     << d_orthogonalization_stop_tol << endl;
  os << "d_orthogonalization_stop_dist = " 
     << d_orthogonalization_stop_dist << endl;
  os << "d_orthogonalization_max_iters = " 
     << d_orthogonalization_max_iters << endl;

  os << "AMR parameters" << endl;
  os << "--------------" << endl;
  os << "d_use_AMR = " << (d_use_AMR ? "true" : "false") << endl;
  os << "d_regrid_interval = " << d_regrid_interval << endl;
  os << "d_tag_buffer_width = " << d_tag_buffer_width << endl;
  os << "d_refinement_cutoff_value = " << d_refinement_cutoff_value << endl;

  os << "PatchData Handles" << endl;
  os << "-----------------" << endl;
//  os << "d_phi_handle = " << d_phi_handle << endl;
// KTC os << "d_phi_time_advance_scr_handles = " << d_phi_time_advance_scr_handles << endl;
//  os << "d_psi_handle = " << d_psi_handle << endl;
// KTC os << "d_psi_time_advance_scr_handles = " << d_psi_time_advance_scr_handles << endl;
  os << "d_grad_phi_upwind_handle = " << d_grad_phi_upwind_handle << endl;
  os << "d_grad_psi_upwind_handle = " << d_grad_psi_upwind_handle << endl;
  os << "d_control_volume_handle = " << d_control_volume_handle << endl;

  os << "Current State" << endl;
  os << "-------------" << endl;
  os << "d_current_time = " << d_current_time << endl;
  os << "d_reinitialization_count = " << d_reinitialization_count << endl;
  os << "d_orthogonalization_count = " << d_orthogonalization_count << endl;
  os << "d_regrid_count = " << d_regrid_count << endl;

  os << "Object Pointers" << endl;
  os << "---------------" << endl;
  os << "d_patch_hierarchy = " << d_patch_hierarchy.getPointer() << endl;
  os << "d_grid_geometry = " << d_grid_geometry.getPointer() << endl;
  os << "d_lsm_patch_strategy = " << d_lsm_patch_strategy << endl; 
  os << "d_lsm_velocity_field_strategy = " << d_lsm_velocity_field_strategy << endl; 
  os << "d_fill_new_level = " << d_fill_new_level.getPointer() << endl; 
  os << "d_fill_bdry_time_advance = " 
     << d_fill_bdry_time_advance.getPointer() << endl; 

  // KTC - add code to print out d_fill_bdry_sched_time_advance

  os << "===================================" << endl << endl;
}


/* computeStableDt() */
template <int DIM> 
LSMLIB_REAL LevelSetFunctionIntegrator<DIM>::computeStableDt()
{
  /*
   * compute maximum stable dt using:
   *
   *   if (user_specified_dt < LSMLIB_REAL_MAX)
   *
   *     dt = user_specified_dt
   *
   *   else
   *
   *     dt = min{physics_dt, advection_dt, normal_vel_dt}
   * 
   * See LevelSetFunctionIntegrator.h for definitions of the four
   * different dt values.  
   * 
   * NOTE: advection_dt and normal_vel_dt are only used if the respective
   *       velocity fields are provided by the 
   *       LevelSetMethodVelocityFieldStrategy.
   */
  LSMLIB_REAL max_advection_dt = LSMLIB_REAL_MAX;
  LSMLIB_REAL max_normal_vel_dt = LSMLIB_REAL_MAX;
  LSMLIB_REAL max_user_specified_dt = LSMLIB_REAL_MAX;

  // allocate scratch space
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<DIM> > level = 
      d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData( d_compute_stable_dt_scratch_variables );
  }
 
  // fill boundary data to for phi/psi to be used for computing
  // velocity field
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: true indicates that physical boundary conditions should
    //       be set.
    d_fill_bdry_sched_compute_stable_dt[ln]->fillData(d_current_time,true);
  }
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
    d_bc_module->imposeBoundaryConditions(
      d_phi_handles[0],
      d_lower_bc_phi[comp], 
      d_upper_bc_phi[comp], 
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      comp);
    if (d_codimension == 2) {
      d_bc_module->imposeBoundaryConditions(
        d_psi_handles[0],
        d_lower_bc_psi[comp], 
        d_upper_bc_psi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
    }
  }

  // loop over PatchHierarchy and compute the maximum stable
  // user-specified dt 
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);

    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name
                  << "::computeStableDt(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      /*
       *  Compute the user-specified dt for the current patch.
       */
      LSMLIB_REAL user_specified_dt_on_patch = 
        d_lsm_patch_strategy->computeStableDtOnPatch(
          *patch,
          this,
          d_lsm_velocity_field_strategy);

      // update max_user_specified_dt
      if ( (max_user_specified_dt > user_specified_dt_on_patch)  &&
           (user_specified_dt_on_patch > 0) ) {
        max_user_specified_dt = user_specified_dt_on_patch;
      }

    }  // end loop over patches in level

  }  // end loop over levels in hierarchy

  max_user_specified_dt = tbox::MPI::minReduction(max_user_specified_dt);

  LSMLIB_REAL max_stable_dt = -1; // temporary value to be set below
  if (max_user_specified_dt < LSMLIB_REAL_MAX) {

    max_stable_dt = max_user_specified_dt;

    if (d_verbose_mode) {
      pout << endl;
      pout << d_object_name << "::computeStableDt():" << endl;
      pout << "  user_specified_dt = " << max_user_specified_dt << endl;
    }

  } else {

    for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

      // compute the velocity field for calculation of 
      // advection_dt and normal_vel_dt
      d_lsm_velocity_field_strategy->computeVelocityField(
        d_current_time, 
        d_phi_handles[0], d_psi_handles[0],
        comp);
  
      /*
       *  If necessary, compute the maximum CFL-based advection dt 
       *  on the current patch
       */
      if (d_lsm_velocity_field_strategy->providesExternalVelocityField()) {
        LSMLIB_REAL max_advection_dt_for_component = 
          LevelSetMethodToolbox<DIM>::computeStableAdvectionDt(
            d_patch_hierarchy,
            d_lsm_velocity_field_strategy->
              getExternalVelocityFieldPatchDataHandle(comp),
            d_control_volume_handle,
            d_cfl_number);

        if ( (max_advection_dt > max_advection_dt_for_component) &&
             (max_advection_dt_for_component > 0) ) {
          max_advection_dt = max_advection_dt_for_component;
        }

      } // end advection velocity dt calculation

      /*
       *  If necessary, compute the maximum CFL-based normal velocity
       *  dt on the current patch
       */
      if (d_lsm_velocity_field_strategy->providesNormalVelocityField()) {
  
        // compute normal velocity dt for phi
        LevelSetMethodToolbox<DIM>::computePlusAndMinusSpatialDerivatives(
          d_patch_hierarchy,
          d_spatial_derivative_type,
          d_spatial_derivative_order,
          d_grad_phi_plus_handle,
          d_grad_phi_minus_handle,
          d_phi_handles[0],
          comp);
  
        LSMLIB_REAL max_phi_normal_vel_dt_for_component = 
          LevelSetMethodToolbox<DIM>::computeStableNormalVelocityDt(
            d_patch_hierarchy,
            d_lsm_velocity_field_strategy->
              getNormalVelocityFieldPatchDataHandle(PHI, comp),
            d_grad_phi_plus_handle,
            d_grad_phi_minus_handle,
            d_control_volume_handle,
            d_cfl_number);
  
        if ( (max_normal_vel_dt > max_phi_normal_vel_dt_for_component) &&
             (max_phi_normal_vel_dt_for_component > 0) ) {
          max_normal_vel_dt = max_phi_normal_vel_dt_for_component;
        }

        if (d_codimension == 2) {
          // compute normal velocity dt for psi
          LevelSetMethodToolbox<DIM>::computePlusAndMinusSpatialDerivatives(
          d_patch_hierarchy,
          d_spatial_derivative_type,
          d_spatial_derivative_order,
          d_grad_psi_plus_handle,
          d_grad_psi_minus_handle,
          d_psi_handles[0],
          comp);
  
          LSMLIB_REAL max_psi_normal_vel_dt_for_component = 
            LevelSetMethodToolbox<DIM>::computeStableNormalVelocityDt(
              d_patch_hierarchy,
              d_lsm_velocity_field_strategy->
                getNormalVelocityFieldPatchDataHandle(PSI, comp),
              d_grad_psi_plus_handle,
              d_grad_psi_minus_handle,
              d_control_volume_handle,
              d_cfl_number);

          if ( (max_normal_vel_dt > max_psi_normal_vel_dt_for_component) &&
               (max_psi_normal_vel_dt_for_component > 0) ) {
            max_normal_vel_dt = max_psi_normal_vel_dt_for_component;
          }
        }
  
      } // end normal velocity dt calculation

    } // end loop over components of level set functions
  
  
    /*
     * compute max_stable_dt
     */
    max_stable_dt = max_advection_dt;

    // take max_normal_vel_dt if it is smaller 
    if (max_normal_vel_dt < max_stable_dt) max_stable_dt = max_normal_vel_dt;
  
    // take physics_dt if it is smaller 
    LSMLIB_REAL physics_dt = d_lsm_velocity_field_strategy->computeStableDt();
    if ( (physics_dt < max_stable_dt) && (physics_dt > 0) ) {
      max_stable_dt = physics_dt;
    }

    /*
     * Take the minimum across all processors.
     * NOTE:  this reduction takes care of two important cases:
     *        (1) max_user_specified_dt is smallest and 
     *        (2) the concrete subclass of LevelSetMethodVelocityFieldStrategy
     *            physics_dt forgets to do a reduction
     */
    max_stable_dt = tbox::MPI::minReduction(max_stable_dt);

    if (d_verbose_mode) {
      pout << endl;
      pout << d_object_name << "::computeStableDt():" << endl;
      pout << "  advection_dt = " << max_advection_dt << endl;
      pout << "  normal_vel_dt = " << max_normal_vel_dt << endl;
      pout << "  physics_dt = " << physics_dt << endl;
    }

  } // end case: user_specified_dt not provided

  // deallocate patch data that was allocated for the time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<DIM> > level 
      = d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData( d_compute_stable_dt_scratch_variables );
  }

  return max_stable_dt;
}


/* advanceLevelSetFunctions() */
template <int DIM> 
bool LevelSetFunctionIntegrator<DIM>::advanceLevelSetFunctions(
  const LSMLIB_REAL dt)
{
  // get number of levels in PatchHierarchy
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  
  // if this is the first time step, synchronize data across processors 
  // NOTE:  normally this is done at the end of the time advance
  if (d_current_time == d_start_time) {
    for ( int ln=0 ; ln < num_levels; ln++ ) {
      // NOTE: true indicates that physical boundary conditions should
      //       be set.
      d_fill_bdry_sched_time_advance[0][ln]
       ->fillData(d_current_time,true);
    }
    for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
      d_bc_module->imposeBoundaryConditions(
        d_phi_handles[0],
        d_lower_bc_phi[comp], 
        d_upper_bc_phi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
      if (d_codimension == 2) {
        d_bc_module->imposeBoundaryConditions(
          d_psi_handles[0],
          d_lower_bc_psi[comp], 
          d_upper_bc_psi[comp], 
          d_spatial_derivative_type,
          d_spatial_derivative_order,
          comp);
      }
    }
  } // end synchronization of data for initial time step

  // allocate scratch space
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<DIM> > level = 
      d_patch_hierarchy->getPatchLevel(ln);
    level->allocatePatchData( d_time_advance_scratch_variables );
  }
 
  // advance level set equation using TVD Runge-Kutta 
  switch(d_tvd_runge_kutta_order) {
    case 1: { // first-order TVD RK (e.g. Forward Euler)
      advanceLevelSetEqnUsingTVDRK1(dt);
      break;
    }
    case 2: { // second-order TVD RK 
      advanceLevelSetEqnUsingTVDRK2(dt);
      break;
    }
    case 3: { // third-order TVD RK 

      advanceLevelSetEqnUsingTVDRK3(dt);
      break;
    }
    default: { // UNSUPPORTED ORDER
      TBOX_ERROR(  d_object_name
                << "::advanceLevelSetFunctions(): " 
                << "Unsupported TVD Runge-Kutta order.  "
                << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
                << endl);
    }
  }

  // increment reinitialization and orthogonalization counters
  d_reinitialization_count++;
  d_orthogonalization_count++;

  // orthogonalize level set functions 
  if ( d_use_orthogonalization &&
       (0 == d_orthogonalization_count % d_orthogonalization_interval) )
  {
    // case: orthogonalization step

    if (d_orthogonalization_evolved_field == PHI) {
      orthogonalizeLevelSetFunctions(PHI);
      d_orthogonalization_evolved_field = PSI;
    } else {
      orthogonalizeLevelSetFunctions(PSI);
      d_orthogonalization_evolved_field = PHI;
    }

    // reset orthogonalization counter 
    d_orthogonalization_count = 0;

    // reset reinitialization counter if necessary
    if ( d_use_reinitialization &&
         (0 == d_reinitialization_count % d_reinitialization_interval) )
    {
      d_reinitialization_count = 0;
    } 

  } else if ( d_use_reinitialization &&
              (0 == d_reinitialization_count % d_reinitialization_interval) )
  {
    // case: reinitialization step, but not an orthogonalization step
    reinitializeLevelSetFunctions(PHI);
    if (d_codimension == 2) {
      reinitializeLevelSetFunctions(PSI);
    }

    // reset reinitialization counter 
    d_reinitialization_count = 0;
  } 

  // determine if patch hierarchy needs to be regridded
  // KTC - DO SOMETHING SMARTER
  // regrid patch hierachy
  bool regrid_needed = false;
  if (0 == d_regrid_count%d_regrid_interval) 
  {
    regrid_needed = true;
    d_regrid_count = 1;
  } else {
    d_regrid_count++;
  } 

  // update current time and integrator step
  d_current_time += dt;
  d_num_integration_steps_taken++;

  // deallocate patch data that was allocated for the time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    Pointer< PatchLevel<DIM> > level 
      = d_patch_hierarchy->getPatchLevel(ln);
    level->deallocatePatchData( d_time_advance_scratch_variables );
  }

  // synchronize data across processors
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: true indicates that physical boundary conditions should
    //       be set.
    d_fill_bdry_sched_time_advance[0][ln]
     ->fillData(d_current_time,true);
  }
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
    d_bc_module->imposeBoundaryConditions(
      d_phi_handles[0],
      d_lower_bc_phi[comp], 
      d_upper_bc_phi[comp], 
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      comp);
    if (d_codimension == 2) {
      d_bc_module->imposeBoundaryConditions(
        d_psi_handles[0],
        d_lower_bc_psi[comp], 
        d_upper_bc_psi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
    }
  }

  // reinitialize level set functions to approximate distance functions
  return regrid_needed;
}


/* preprocessInitializeVelocityField() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::preprocessInitializeVelocityField(
  int& phi_handle,
  int& psi_handle,
  const Pointer< PatchHierarchy<DIM> > hierarchy,
  const int level_number)
{
  Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(level_number);

  // fill scratch data
  Pointer< RefineSchedule<DIM> > sched =
    d_fill_new_level->createSchedule(level, level_number-1, hierarchy, this);
  sched->fillData(d_current_time,true);

  phi_handle = d_phi_handles[0]; 
  psi_handle = d_psi_handles[0]; 
}


/* postprocessInitializeVelocityField() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::postprocessInitializeVelocityField(
  const Pointer< PatchHierarchy<DIM> > hierarchy,
  const int level_number)
{
}


/* putToDatabase() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::putToDatabase(Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   // write out version information
   db->putInteger("LEVEL_SET_METHOD_INTEGRATOR_VERSION",
                   LEVEL_SET_METHOD_INTEGRATOR_VERSION);

  /*
   * Write user-specified parameters to database
   */

  // read in level set parameters
  db->putDouble("d_start_time", d_start_time);
  db->putDouble("d_end_time", d_end_time);
  db->putDouble("d_cfl_number", d_cfl_number);

  db->putInteger("d_spatial_derivative_type", d_spatial_derivative_type);
  db->putInteger("d_spatial_derivative_order", d_spatial_derivative_order);
  db->putInteger("d_tvd_runge_kutta_order", d_tvd_runge_kutta_order);

  db->putInteger("d_reinitialization_interval", d_reinitialization_interval);
  db->putDouble("d_reinitialization_stop_tol", d_reinitialization_stop_tol);
  db->putDouble("d_reinitialization_stop_dist", d_reinitialization_stop_dist);
  db->putInteger("d_reinitialization_max_iters", d_reinitialization_max_iters);
  db->putInteger("d_orthogonalization_interval", d_orthogonalization_interval);
  db->putDouble("d_orthogonalization_stop_tol", d_orthogonalization_stop_tol);
  db->putDouble("d_orthogonalization_stop_dist", 
                  d_orthogonalization_stop_dist);
  db->putInteger("d_orthogonalization_max_iters", 
                  d_orthogonalization_max_iters);

  db->putBool("d_use_AMR", d_use_AMR); 
  db->putInteger("d_regrid_interval", d_regrid_interval);
  db->putInteger("d_tag_buffer_width", d_tag_buffer_width);
  db->putDouble("d_refinement_cutoff_value", d_refinement_cutoff_value);

  db->putBool("d_verbose_mode", d_verbose_mode); 

  /*
   * Write state parameters to database
   */
  db->putDouble("d_current_time", d_current_time);
  db->putInteger("d_num_integration_steps_taken", 
    d_num_integration_steps_taken);
  db->putInteger("d_reinitialization_count", d_reinitialization_count);
  db->putInteger("d_orthogonalization_count", d_orthogonalization_count);
  db->putInteger("d_regrid_count", d_regrid_count);

  db->putBool("d_use_reinitialization", d_use_reinitialization);
  db->putBool("d_use_reinitialization_stop_tol", 
    d_use_reinitialization_stop_tol);
  db->putBool("d_use_reinitialization_stop_dist", 
    d_use_reinitialization_stop_dist);
  db->putBool("d_use_reinitialization_max_iters", 
    d_use_reinitialization_max_iters);
  db->putBool("d_use_orthogonalization", d_use_orthogonalization);
  db->putBool("d_use_orthogonalization_stop_tol", 
    d_use_orthogonalization_stop_tol);
  db->putBool("d_use_orthogonalization_stop_dist", 
    d_use_orthogonalization_stop_dist);
  db->putBool("d_use_orthogonalization_max_iters", 
    d_use_orthogonalization_max_iters);

}


/* initializeLevelData() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::initializeLevelData (
  const Pointer< BasePatchHierarchy<DIM> > hierarchy ,
  const int level_number ,
  const double init_data_time ,
  const bool can_be_refined ,
  const bool initial_time ,
  const Pointer< BasePatchLevel<DIM> > old_level,
  const bool allocate_data)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!hierarchy.isNull());
   assert( (level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()) );
   if ( !(old_level.isNull()) ) {
      assert( level_number == old_level->getLevelNumber() );
   }
   assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

  Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(level_number);

  if (allocate_data) {
    level->allocatePatchData(d_solution_variables);
    level->allocatePatchData(d_persistent_variables);
  } else {
    level->setTime( init_data_time, d_solution_variables);
    level->setTime( init_data_time, d_persistent_variables);
  }

  // create schedules for filling new level and fill new level with data
  if ((level_number > 0) || !old_level.isNull()) {

    Pointer< RefineSchedule<DIM> > sched =
      d_fill_new_level->createSchedule(
        level, old_level, level_number-1,
        hierarchy, this);

    sched->fillData(init_data_time);

  }

  /*
   * Manually initialize data on all patches in the level if this
   * is the initial time.
   */ 
  if (initial_time)
  {
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) {
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name 
                  << "::initializeLevelData(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl );
      }

      d_lsm_patch_strategy->initializeLevelSetFunctionsOnPatch(
        *patch, init_data_time, d_phi_handles[0], d_psi_handles[0]);

    } // end loop over patches
  } // end (initial_time)

  // deallocate data on old PatchLevel if it exists
  if (!old_level.isNull()) {
    old_level->deallocatePatchData(d_solution_variables);
    old_level->deallocatePatchData(d_persistent_variables);
  }

}


/* setPhysicalBoundaryConditions() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::setPhysicalBoundaryConditions(
  Patch<DIM>& patch,
  const double fill_time,
  const IntVector<DIM>& ghost_width_to_fill)
{
  // invoke LevelSetMethodPatchStrategy::setPhysicalBoundaryConditions()
  d_lsm_patch_strategy
    ->setLevelSetFunctionBoundaryConditions(patch, fill_time, 
                                            d_phi_handles[0], d_psi_handles[0], 
                                            ghost_width_to_fill);
}


/* setBoundaryConditions() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::setBoundaryConditions(
  const IntVector<DIM>& lower_bc,
  const IntVector<DIM>& upper_bc,
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int component)
{

  // check that anti-periodic boundary conditions are consistent
  for (int dim = 0; dim < DIM; dim++) {
    if ( ((lower_bc[dim] == BoundaryConditionModule<DIM>::ANTI_PERIODIC) &&
          (upper_bc[dim] != BoundaryConditionModule<DIM>::ANTI_PERIODIC)) ||
         ((lower_bc[dim] != BoundaryConditionModule<DIM>::ANTI_PERIODIC) &&
          (upper_bc[dim] == BoundaryConditionModule<DIM>::ANTI_PERIODIC)) ) {
      TBOX_ERROR(  d_object_name 
                << "::setBoundaryConditions(): "
                << "ANTI_PERIODIC boundary conditions inconsistently set."
                << endl );
    }
  }

  if (level_set_fcn == LSMLIB::PHI) {
    if (component >= 0) {
      d_lower_bc_phi[component] = lower_bc;
      d_upper_bc_phi[component] = upper_bc;
    } else {
      for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
        d_lower_bc_phi[comp] = lower_bc;
        d_upper_bc_phi[comp] = upper_bc;
      } 
    } 
  }
  else if ( (d_codimension == 2) && (level_set_fcn == LSMLIB::PSI) ) {
    if (component >= 0) {
      d_lower_bc_psi[component] = lower_bc;
      d_upper_bc_psi[component] = upper_bc;
    } else {
      for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
        d_lower_bc_psi[comp] = lower_bc;
        d_upper_bc_psi[comp] = upper_bc;
      } 
    } 
  }
}


/* applyGradientDetector() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::applyGradientDetector(
  const Pointer< BasePatchHierarchy<DIM> > hierarchy,
  const int level_number,
  const double error_data_time,
  const int tag_index,
  const bool initial_time,
  const bool uses_richardson_extrapolation_too)
{
  Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(level_number);
  typename PatchLevel<DIM>::Iterator pi;
  for (pi.initialize(level); pi; pi++) { // loop over patches
    const int pn = *pi;
    Pointer< Patch<DIM> > patch = level->getPatch(pn);
    if ( patch.isNull() ) {
      TBOX_ERROR(  d_object_name 
                << "::applyGradientDetector(): "
                << "Cannot find patch. Null patch pointer."
                << endl );
    }

    Pointer< CellData<DIM,LSMLIB_REAL> > phi_data = 
      patch->getPatchData( d_phi_handles[0] );
    Pointer< CellData<DIM,LSMLIB_REAL> > psi_data;
    if (d_codimension == 2) {
      psi_data = patch->getPatchData( d_psi_handles[0] );
    }
    Pointer< CellData<DIM,int> > tag_data = 
      patch->getPatchData( tag_index );

    Box<DIM> phi_ghostbox = phi_data->getGhostBox();
    const IntVector<DIM> phi_gb_lower = phi_ghostbox.lower();
    const IntVector<DIM> phi_gb_upper = phi_ghostbox.upper();

    Box<DIM> psi_ghostbox = psi_data->getGhostBox();
    const IntVector<DIM> psi_gb_lower = psi_ghostbox.lower();
    const IntVector<DIM> psi_gb_upper = psi_ghostbox.upper();

    Box<DIM> tag_box = tag_data->getBox();
    const IntVector<DIM> tag_box_lower = tag_box.lower();
    const IntVector<DIM> tag_box_upper = tag_box.upper();

    Box<DIM> tag_ghostbox = tag_data->getGhostBox();
    const IntVector<DIM> tag_gb_lower = tag_ghostbox.lower();
    const IntVector<DIM> tag_gb_upper = tag_ghostbox.upper();

/* KTC - FIX ME
#if DIM==2
    setrefinementtags_(
      tag_data->getPointer(),
      level_set_data->getPointer(),
      &tag_ghostbox_lower[0],
      &tag_ghostbox_upper[0],
      &tag_ghostbox_lower[1],
      &tag_ghostbox_upper[1],
      &level_set_ghostbox_lower[0],
      &level_set_ghostbox_upper[0],
      &level_set_ghostbox_lower[1],
      &level_set_ghostbox_upper[1],
      &tag_box_lower[0],
      &tag_box_upper[0],
      &tag_box_lower[1],
      &tag_box_upper[1],
      &d_refinement_cutoff_value);
#endif

#if DIM==3
#endif
*/

  } // end loop over patches

}


/* resetHierarchyConfiguration() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::resetHierarchyConfiguration(
  Pointer< BasePatchHierarchy<DIM> > hierarchy ,
  int coarsest_level ,
  int finest_level )
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!hierarchy.isNull());
  assert( (coarsest_level >= 0)
          && (coarsest_level <= finest_level)
          && (finest_level <= hierarchy->getFinestLevelNumber()) );
  for (int ln = 0; ln <= finest_level; ln++) {
     assert(!(hierarchy->getPatchLevel(ln)).isNull());
  }
#endif

  int num_levels = hierarchy->getNumberLevels();

  // reset communications schedules used to fill boundary data 
  // during time advance
  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_fill_bdry_sched_time_advance[k].resizeArray(num_levels);

    for (int ln = coarsest_level; ln <= finest_level; ln++) {
      Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);
  
      // reset data transfer configuration for boundary filling 
      // before time advance
      d_fill_bdry_sched_time_advance[k][ln] =
        d_fill_bdry_time_advance[k]->createSchedule(level,
                                                    ln-1,
                                                    hierarchy,
                                                    this);
  
    } // end loop over levels
  } // end loop over TVD Runge-Kutta stages

  // reset communications schedules used to fill boundary data when
  // computing the stable dt 
  d_fill_bdry_sched_compute_stable_dt.resizeArray(num_levels);

  for (int ln = coarsest_level; ln <= finest_level; ln++) {
    Pointer< PatchLevel<DIM> > level = hierarchy->getPatchLevel(ln);

    // reset data transfer configuration for boundary filling 
    // before time advance
    d_fill_bdry_sched_compute_stable_dt[ln] =
      d_fill_bdry_compute_stable_dt->createSchedule(level,
                                                    ln-1,
                                                    hierarchy,
                                                    this);
  } // end loop over levels

  // reset hierarchy configuration for reinitialization and orthogonalization
  // algorithms
  d_phi_reinitialization_alg->resetHierarchyConfiguration(
    hierarchy, coarsest_level, finest_level);
  if (d_codimension == 2) {
    d_psi_reinitialization_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);
    d_orthogonalization_alg->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level);
  }

  // recompute control volumes
  LevelSetMethodToolbox<DIM>::computeControlVolumes(
    hierarchy, d_control_volume_handle);

  // create anti-periodic boundary condition plan
  d_bc_module->resetHierarchyConfiguration(
    d_patch_hierarchy,
    coarsest_level,
    finest_level,
    d_level_set_ghostcell_width);

}


/* advanceLevelSetEqnUsingTVDRK1() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::advanceLevelSetEqnUsingTVDRK1(
  const LSMLIB_REAL dt)
{
  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time, 
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through TVD-RK1 step 
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK1Step(
      d_patch_hierarchy,
      d_phi_handles[0], 
      d_phi_handles[rk_stage], 
      d_rhs_phi_handle, dt,
      comp, comp, 0); // components of PatchData to use in TVD-RK1 step

    if (d_codimension == 2) {

      // advance psi through TVD-RK1 step 
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK1Step(
        d_patch_hierarchy,
        d_psi_handles[0], 
        d_psi_handles[rk_stage], 
        d_rhs_psi_handle, dt,
        comp, comp, 0); // components of PatchData to use in TVD-RK1 step

    } // end codimension-two case

  } // end loop over vector level set function
}


/* advanceLevelSetEqnUsingTVDRK2() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::advanceLevelSetEqnUsingTVDRK2(
  const LSMLIB_REAL dt)
{
  // { begin Stage 1

  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time, 
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through the first stage of TVD-RK2 
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK2Stage1(
      d_patch_hierarchy,
      d_phi_handles[rk_stage+1],
      d_phi_handles[rk_stage],
      d_rhs_phi_handle, dt,
      comp, comp, 0); // components of PatchData to use in first 
                      // stage TVD-RK2 step

    if (d_codimension == 2) {

      // advance psi through the first stage of TVD-RK2 
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK2Stage1(
        d_patch_hierarchy,
        d_psi_handles[rk_stage+1],
        d_psi_handles[rk_stage],
        d_rhs_psi_handle, dt,
        comp, comp, 0); // components of PatchData to use in first 
                        // stage TVD-RK2 step
    }
  } // end loop over vector level set function

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK2 stage counter
  rk_stage = 1;

  // fill scratch space for second stage of time advance
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: true indicates that physical boundary conditions should
    //       be set.
    d_fill_bdry_sched_time_advance[rk_stage][ln]
     ->fillData(d_current_time,true);
  }
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
    d_bc_module->imposeBoundaryConditions(
      d_phi_handles[rk_stage],
      d_lower_bc_phi[comp], 
      d_upper_bc_phi[comp], 
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      comp);
    if (d_codimension == 2) {
      d_bc_module->imposeBoundaryConditions(
        d_psi_handles[rk_stage],
        d_lower_bc_psi[comp], 
        d_upper_bc_psi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
    }
  }

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time+dt, 
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through the second stage of TVD-RK2 
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK2Stage2(
      d_patch_hierarchy,
      d_phi_handles[0],
      d_phi_handles[rk_stage],
      d_phi_handles[0],
      d_rhs_phi_handle, dt,
      comp, comp, comp, 0); // components of PatchData to use in final 
                            // stage of TVD-RK2 step

    if (d_codimension == 2) {

      // advance psi through the second stage of TVD-RK2 
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK2Stage2(
        d_patch_hierarchy,
        d_psi_handles[0],
        d_psi_handles[rk_stage],
        d_psi_handles[0],
        d_rhs_psi_handle, dt,
        comp, comp, comp, 0); // components of PatchData to use in final 
                              // stage of TVD-RK2 step
    }
  } // end loop over components of vector level set function

  // } end Stage 2
}


/* advanceLevelSetEqnUsingTVDRK3() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::advanceLevelSetEqnUsingTVDRK3(
  const LSMLIB_REAL dt)
{
  // { begin Stage 1

  // initialize counter for current stage of TVD RK step
  // NOTE: the rk_stage begins at 0 for convenience
  int rk_stage = 0;

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time,
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through the first stage of TVD-RK3
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK3Stage1(
      d_patch_hierarchy,
      d_phi_handles[rk_stage+1],
      d_phi_handles[rk_stage],
      d_rhs_phi_handle, dt,
      comp, comp, 0); // components of PatchData to use in first stage 
                      // of TVD-RK3 step

    if (d_codimension == 2) {
  
      // advance psi through the first stage of TVD-RK3
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK3Stage1(
        d_patch_hierarchy,
        d_psi_handles[rk_stage+1],
        d_psi_handles[rk_stage],
        d_rhs_psi_handle, dt,
        comp, comp, 0); // components of PatchData to use in first stage 
                        // of TVD-RK3 step
    }
  } // end loop over vector level set function

  // } end Stage 1


  // { begin Stage 2

  // advance TVD RK3 stage counter
  rk_stage = 1;

  // fill scratch space for second stage of time advance
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: true indicates that physical boundary conditions should
    //       be set.
    d_fill_bdry_sched_time_advance[rk_stage][ln]
     ->fillData(d_current_time,true);
  }
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
    d_bc_module->imposeBoundaryConditions(
      d_phi_handles[rk_stage],
      d_lower_bc_phi[comp], 
      d_upper_bc_phi[comp], 
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      comp);
    if (d_codimension == 2) {
      d_bc_module->imposeBoundaryConditions(
        d_psi_handles[rk_stage],
        d_lower_bc_psi[comp], 
        d_upper_bc_psi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
    }
  }

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time+dt,
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through the second stage of TVD-RK3
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK3Stage2(
      d_patch_hierarchy,
      d_phi_handles[rk_stage+1],
      d_phi_handles[rk_stage],
      d_phi_handles[rk_stage-1],
      d_rhs_phi_handle, dt,
      comp, comp, comp, 0); // components of PatchData to use in second 
                            // stage of TVD-RK3 step

    if (d_codimension == 2) {

      // advance psi through the second stage of TVD-RK3
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK3Stage2(
        d_patch_hierarchy,
        d_psi_handles[rk_stage+1],
        d_psi_handles[rk_stage],
        d_psi_handles[rk_stage-1],
        d_rhs_psi_handle, dt,
        comp, comp, comp, 0); // components of PatchData to use in second 
                              // stage of TVD-RK3 step
    }
  } // end loop over vector level set function

  // } end Stage 2


  // { begin Stage 3

  // advance TVD RK3 stage counter
  rk_stage = 2;

  // fill scratch space for second stage of time advance
  for ( int ln=0 ; ln < num_levels; ln++ ) {
    // NOTE: true indicates that physical boundary conditions should
    //       be set.
    d_fill_bdry_sched_time_advance[rk_stage][ln]
     ->fillData(d_current_time,true);
  }
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {
    d_bc_module->imposeBoundaryConditions(
      d_phi_handles[rk_stage],
      d_lower_bc_phi[comp], 
      d_upper_bc_phi[comp], 
      d_spatial_derivative_type,
      d_spatial_derivative_order,
      comp);
    if (d_codimension == 2) {
      d_bc_module->imposeBoundaryConditions(
        d_psi_handles[rk_stage],
        d_lower_bc_psi[comp], 
        d_upper_bc_psi[comp], 
        d_spatial_derivative_type,
        d_spatial_derivative_order,
        comp);
    }
  }

  // loop over components of vector level set function
  for (int comp = 0; comp < d_num_level_set_fcn_components; comp++) {

    // compute velocity field for current stage
    d_lsm_velocity_field_strategy->computeVelocityField(
      d_current_time+0.5*dt,
      d_phi_handles[rk_stage],
      d_psi_handles[rk_stage],
      comp);

    // advance phi through the second stage of TVD-RK3
    computeLevelSetEquationRHS(PHI,d_phi_handles[rk_stage],
                               comp);
    LevelSetMethodToolbox<DIM>::TVDRK3Stage3(
      d_patch_hierarchy,
      d_phi_handles[0],
      d_phi_handles[rk_stage],
      d_phi_handles[0],
      d_rhs_phi_handle, dt,
      comp, comp, comp, 0); // components of PatchData to use in final 
                            // stage of TVD-RK3 step

    if (d_codimension == 2) {
  
      // advance psi through the second stage of TVD-RK3
      computeLevelSetEquationRHS(PSI,d_psi_handles[rk_stage],
                                 comp);
      LevelSetMethodToolbox<DIM>::TVDRK3Stage3(
        d_patch_hierarchy,
        d_psi_handles[0],
        d_psi_handles[rk_stage],
        d_psi_handles[0],
        d_rhs_psi_handle, dt,
        comp, comp, comp, 0); // components of PatchData to use in final 
                              // stage of TVD-RK3 step
    }
  } // end loop over vector level set function
  
  // } end Stage 3
}


/* computeLevelSetEquationRHS() first zeros out the RHS and then
 * calls addAdvectionTermToLevelSetEquationRHS() and 
 * addNormalVelocityTermToLevelSetEquationRHS() as appropriate.
 */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::computeLevelSetEquationRHS(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int phi_handle,
  const int component)
{
  int rhs_handle;
  if (level_set_fcn == PHI) {
    rhs_handle = d_rhs_phi_handle;
  } else {
    rhs_handle = d_rhs_psi_handle;
  } 

  // loop over PatchHierarchy and zero out the RHS for level set 
  // equation by calling Fortran routines
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name 
                  << "::computeLevelSetEquationRHS(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,LSMLIB_REAL> > rhs_data =
        patch->getPatchData( rhs_handle );
  
      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      LSMLIB_REAL* rhs = rhs_data->getPointer();

      // zero out level set equation RHS
      if (DIM == 3) {

        LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2]);

      } else if (DIM == 2) {

        LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1]);

      } else if (DIM == 1) {

        LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  d_object_name 
                  << "::computeLevelSetEquationRHS(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      } // end switch over dimension (DIM) of calculation

    } // end loop over patches in level
  } // end loop over levels in hierarchy

  // invoke addAdvectionTermToLevelSetEquationRHS() if necessary
  if (d_lsm_velocity_field_strategy->providesExternalVelocityField()) {
    addAdvectionTermToLevelSetEquationRHS(level_set_fcn, phi_handle,
                                          component);
  }

  // invoke addNormalVelocityTermToLevelSetEquationRHS() if necessary
  if (d_lsm_velocity_field_strategy->providesNormalVelocityField()) {
    addNormalVelocityTermToLevelSetEquationRHS(level_set_fcn, phi_handle,
                                               component);
  }

}


/* addAdvectionTermToLevelSetEquationRHS() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::addAdvectionTermToLevelSetEquationRHS(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int phi_handle,
  const int component)
{
  int grad_phi_upwind_handle;
  int rhs_handle;
  if (level_set_fcn == PHI) {
    grad_phi_upwind_handle = d_grad_phi_upwind_handle;
    rhs_handle = d_rhs_phi_handle;
  } else {
    grad_phi_upwind_handle = d_grad_psi_upwind_handle;
    rhs_handle = d_rhs_psi_handle;
  } 

  int velocity_handle = d_lsm_velocity_field_strategy->
    getExternalVelocityFieldPatchDataHandle(component);

  // compute upwind spatial derivatives 
  LevelSetMethodToolbox<DIM>::computeUpwindSpatialDerivatives(
    d_patch_hierarchy,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    grad_phi_upwind_handle,
    phi_handle,
    velocity_handle,
    component); 

  // loop over PatchHierarchy and add contribution of advection term 
  // to level set equation RHS by calling Fortran subroutines
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name 
                  << "::addAdvectionTermToLevelSetEquationRHS(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,LSMLIB_REAL> > rhs_data =
        patch->getPatchData( rhs_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > grad_phi_upwind_data =
        patch->getPatchData( grad_phi_upwind_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > velocity_data =
        patch->getPatchData( velocity_handle );

      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      Box<DIM> grad_phi_upwind_ghostbox = grad_phi_upwind_data->getGhostBox();
      const IntVector<DIM> grad_phi_upwind_ghostbox_lower = 
        grad_phi_upwind_ghostbox.lower();
      const IntVector<DIM> grad_phi_upwind_ghostbox_upper = 
        grad_phi_upwind_ghostbox.upper();

      Box<DIM> vel_ghostbox = velocity_data->getGhostBox();
      const IntVector<DIM> vel_ghostbox_lower = vel_ghostbox.lower();
      const IntVector<DIM> vel_ghostbox_upper = vel_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = rhs_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      LSMLIB_REAL* rhs = rhs_data->getPointer();
      LSMLIB_REAL* grad_phi_upwind[LSM_DIM_MAX];
      LSMLIB_REAL* vel[LSM_DIM_MAX];
      for (int dim = 0; dim < DIM; dim++) {
        grad_phi_upwind[dim] = grad_phi_upwind_data->getPointer(dim);
        vel[dim] = velocity_data->getPointer(dim);
      }

      if (DIM == 3) {

        LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          grad_phi_upwind[0], grad_phi_upwind[1], grad_phi_upwind[2],
          &grad_phi_upwind_ghostbox_lower[0],
          &grad_phi_upwind_ghostbox_upper[0],
          &grad_phi_upwind_ghostbox_lower[1],
          &grad_phi_upwind_ghostbox_upper[1],
          &grad_phi_upwind_ghostbox_lower[2],
          &grad_phi_upwind_ghostbox_upper[2],
          vel[0], vel[1], vel[2],
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2]);

      } else if (DIM == 2) {

        LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          grad_phi_upwind[0], grad_phi_upwind[1],
          &grad_phi_upwind_ghostbox_lower[0],
          &grad_phi_upwind_ghostbox_upper[0],
          &grad_phi_upwind_ghostbox_lower[1],
          &grad_phi_upwind_ghostbox_upper[1],
          vel[0], vel[1],
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1]);

      } else if (DIM == 1) {

        LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          grad_phi_upwind[0], 
          &grad_phi_upwind_ghostbox_lower[0],
          &grad_phi_upwind_ghostbox_upper[0],
          vel[0], 
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  d_object_name 
                  << "::addAdvectionTermToLevelSetEquationRHS(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      } // end switch over dimension (DIM) of calculation

    } // end loop over patches in level
  } // end loop over levels in hierarchy
}


/* addNormalVelocityTermToLevelSetEquationRHS() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::addNormalVelocityTermToLevelSetEquationRHS(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int phi_handle,
  const int component)
{
  int grad_phi_plus_handle;
  int grad_phi_minus_handle;
  int rhs_handle;
  if (level_set_fcn == PHI) {
    grad_phi_plus_handle = d_grad_phi_plus_handle;
    grad_phi_minus_handle = d_grad_phi_minus_handle;
    rhs_handle = d_rhs_phi_handle;
  } else {
    grad_phi_plus_handle = d_grad_psi_plus_handle;
    grad_phi_minus_handle = d_grad_psi_minus_handle;
    rhs_handle = d_rhs_psi_handle;
  } 

  int normal_velocity_handle = d_lsm_velocity_field_strategy->
    getNormalVelocityFieldPatchDataHandle(level_set_fcn, component);

  // compute plus and minus spatial derivatives 
  LevelSetMethodToolbox<DIM>::computePlusAndMinusSpatialDerivatives(
    d_patch_hierarchy,
    d_spatial_derivative_type,
    d_spatial_derivative_order,
    grad_phi_plus_handle,
    grad_phi_minus_handle,
    phi_handle,
    component); 

  // loop over PatchHierarchy and add contribution of normal velocity 
  // term to level set equation RHS by calling Fortran subroutines
  const int num_levels = d_patch_hierarchy->getNumberLevels();
  for ( int ln=0 ; ln < num_levels; ln++ ) {

    Pointer< PatchLevel<DIM> > level = d_patch_hierarchy->getPatchLevel(ln);
    
    typename PatchLevel<DIM>::Iterator pi;
    for (pi.initialize(level); pi; pi++) { // loop over patches
      const int pn = *pi;
      Pointer< Patch<DIM> > patch = level->getPatch(pn);
      if ( patch.isNull() ) {
        TBOX_ERROR(  d_object_name 
                  << "::addNormalVelocityTermToLevelSetEquationRHS(): "
                  << "Cannot find patch. Null patch pointer."
                  << endl);
      }

      // get pointers to data and index space ranges
      Pointer< CellData<DIM,LSMLIB_REAL> > rhs_data =
        patch->getPatchData( rhs_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > normal_velocity_data =
        patch->getPatchData( normal_velocity_handle );

      Pointer< CellData<DIM,LSMLIB_REAL> > grad_phi_plus_data =
        patch->getPatchData( grad_phi_plus_handle );
      Pointer< CellData<DIM,LSMLIB_REAL> > grad_phi_minus_data =
        patch->getPatchData( grad_phi_minus_handle );
  
      Box<DIM> rhs_ghostbox = rhs_data->getGhostBox();
      const IntVector<DIM> rhs_ghostbox_lower = rhs_ghostbox.lower();
      const IntVector<DIM> rhs_ghostbox_upper = rhs_ghostbox.upper();

      Box<DIM> grad_phi_plus_ghostbox = grad_phi_plus_data->getGhostBox();
      const IntVector<DIM> grad_phi_plus_ghostbox_lower = 
        grad_phi_plus_ghostbox.lower();
      const IntVector<DIM> grad_phi_plus_ghostbox_upper = 
        grad_phi_plus_ghostbox.upper();
      Box<DIM> grad_phi_minus_ghostbox = grad_phi_minus_data->getGhostBox();
      const IntVector<DIM> grad_phi_minus_ghostbox_lower = 
        grad_phi_minus_ghostbox.lower();
      const IntVector<DIM> grad_phi_minus_ghostbox_upper = 
        grad_phi_minus_ghostbox.upper();

      Box<DIM> vel_ghostbox = normal_velocity_data->getGhostBox();
      const IntVector<DIM> vel_ghostbox_lower = vel_ghostbox.lower();
      const IntVector<DIM> vel_ghostbox_upper = vel_ghostbox.upper();

      // fill box
      Box<DIM> fillbox = rhs_data->getBox();
      const IntVector<DIM> fillbox_lower = fillbox.lower();
      const IntVector<DIM> fillbox_upper = fillbox.upper();

      LSMLIB_REAL* rhs = rhs_data->getPointer();
      LSMLIB_REAL* grad_phi_plus[LSM_DIM_MAX];
      LSMLIB_REAL* grad_phi_minus[LSM_DIM_MAX];
      LSMLIB_REAL* vel = normal_velocity_data->getPointer();
      for (int dim = 0; dim < DIM; dim++) {
        grad_phi_plus[dim] = grad_phi_plus_data->getPointer(dim);
        grad_phi_minus[dim] = grad_phi_minus_data->getPointer(dim);
      }

      if (DIM == 3) {

        LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          &rhs_ghostbox_lower[2],
          &rhs_ghostbox_upper[2],
          grad_phi_plus[0], grad_phi_plus[1], grad_phi_plus[2],
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          &grad_phi_plus_ghostbox_lower[2],
          &grad_phi_plus_ghostbox_upper[2],
          grad_phi_minus[0], grad_phi_minus[1], grad_phi_minus[2],
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          &grad_phi_minus_ghostbox_lower[2],
          &grad_phi_minus_ghostbox_upper[2],
          vel,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &vel_ghostbox_lower[2],
          &vel_ghostbox_upper[2],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1],
          &fillbox_lower[2],
          &fillbox_upper[2]);

      } else if (DIM == 2) {

        LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          &rhs_ghostbox_lower[1],
          &rhs_ghostbox_upper[1],
          grad_phi_plus[0], grad_phi_plus[1],
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          &grad_phi_plus_ghostbox_lower[1],
          &grad_phi_plus_ghostbox_upper[1],
          grad_phi_minus[0], grad_phi_minus[1],
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          &grad_phi_minus_ghostbox_lower[1],
          &grad_phi_minus_ghostbox_upper[1],
          vel,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &vel_ghostbox_lower[1],
          &vel_ghostbox_upper[1],
          &fillbox_lower[0],
          &fillbox_upper[0],
          &fillbox_lower[1],
          &fillbox_upper[1]);

      } else if (DIM == 1) {

        LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS(
          rhs,
          &rhs_ghostbox_lower[0],
          &rhs_ghostbox_upper[0],
          grad_phi_plus[0], 
          &grad_phi_plus_ghostbox_lower[0],
          &grad_phi_plus_ghostbox_upper[0],
          grad_phi_minus[0], 
          &grad_phi_minus_ghostbox_lower[0],
          &grad_phi_minus_ghostbox_upper[0],
          vel,
          &vel_ghostbox_lower[0],
          &vel_ghostbox_upper[0],
          &fillbox_lower[0],
          &fillbox_upper[0]);

      } else {  // Unsupported dimension
        TBOX_ERROR(  d_object_name 
                  << "::addNormalVelocityTermToLevelSetEquationRHS(): "
                  << "Invalid value of DIM.  "
                  << "Only DIM = 1, 2, and 3 are supported."
                  << endl);
      } // end switch over dimension (DIM) of calculation

    } // end loop over patches in level
  } // end loop over levels in hierarchy
}


/* reinitializeLevelSetFunctions() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::reinitializeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int max_iterations)
{
  if (level_set_fcn == PHI) {
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      d_phi_reinitialization_alg->
        reinitializeLevelSetFunctionForSingleComponent(
          comp,
          max_iterations,
          d_lower_bc_phi[comp],
          d_upper_bc_phi[comp]);
    }
  } else {
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      d_psi_reinitialization_alg->
        reinitializeLevelSetFunctionForSingleComponent(
          comp,
          max_iterations,
          d_lower_bc_psi[comp],
          d_upper_bc_psi[comp]);
    }
  }
}


/* orthogonalizeLevelSetFunctions() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::orthogonalizeLevelSetFunctions(
  const LEVEL_SET_FCN_TYPE level_set_fcn,
  const int max_reinit_iterations,
  const int max_ortho_iterations)
{
  if (level_set_fcn == PHI) {
    reinitializeLevelSetFunctions(PSI, max_reinit_iterations);
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      d_orthogonalization_alg->
        orthogonalizeLevelSetFunctionForSingleComponent(PHI, 
          comp,
          max_ortho_iterations,
          d_lower_bc_psi[comp],
          d_upper_bc_psi[comp],
          d_lower_bc_phi[comp],
          d_upper_bc_phi[comp]);
    }
    reinitializeLevelSetFunctions(PHI, max_reinit_iterations);
  } else {
    reinitializeLevelSetFunctions(PHI, max_reinit_iterations);
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      d_orthogonalization_alg->
        orthogonalizeLevelSetFunctionForSingleComponent(PSI, 
          comp,
          max_ortho_iterations,
          d_lower_bc_phi[comp],
          d_upper_bc_phi[comp],
          d_lower_bc_psi[comp],
          d_upper_bc_psi[comp]);
    }
    reinitializeLevelSetFunctions(PSI, max_reinit_iterations);
  }
}


/* initializeVariables() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::initializeVariables()
{
  // setup ghost cell widths
  int scratch_ghostcell_width = -1; // bogus value set in switch statement
  switch (d_spatial_derivative_type) {
    case ENO: {
      scratch_ghostcell_width = d_spatial_derivative_order;
      break;
    }
    case WENO: {
      scratch_ghostcell_width = d_spatial_derivative_order/2 + 1;
      break;
    }
    default:
      TBOX_ERROR(d_object_name 
              << "::initializeVariables(): "
              << "Unsupported spatial derivative type.  "
              << "Only ENO and WENO derivatives are supported."
              << endl );
  }

  d_level_set_ghostcell_width = IntVector<DIM>(scratch_ghostcell_width);
  IntVector<DIM> zero_ghostcell_width(0);

  // get pointer to VariableDatabase
  VariableDatabase<DIM> *var_db = VariableDatabase<DIM>::getDatabase();

  // get contexts used by level set method algorithm
  Pointer<VariableContext> current_context = var_db->getContext("CURRENT");
  Pointer<VariableContext> scratch_context = var_db->getContext("SCRATCH");
  Pointer<VariableContext> upwind_context = var_db->getContext("UPWIND");
  Pointer<VariableContext> plus_context = var_db->getContext("PLUS_DERIVATIVE");
  Pointer<VariableContext> minus_context = 
    var_db->getContext("MINUS_DERIVATIVE");

  // initialize LevelSetFunctionIntegrator component selectors
  d_solution_variables.clrAllFlags();
  d_compute_stable_dt_scratch_variables.clrAllFlags();
  d_reinitialization_scratch_variables.clrAllFlags();
  d_orthogonalization_scratch_variables.clrAllFlags();
  d_time_advance_scratch_variables.clrAllFlags();
  d_persistent_variables.clrAllFlags();


  /* 
   * Initialize phi variables
   */

  // create variables
  Pointer< CellVariable<DIM,LSMLIB_REAL> > phi_variable;
  if (var_db->checkVariableExists("phi (LSMLIB)")) {
    phi_variable = var_db->getVariable("phi (LSMLIB)");
  } else {
    phi_variable = new CellVariable<DIM,LSMLIB_REAL>(
      "phi (LSMLIB)", d_num_level_set_fcn_components); 
  }
  Pointer< CellVariable<DIM,LSMLIB_REAL> > grad_phi_variable;
  if (var_db->checkVariableExists("grad phi (LSMLIB)")) {
    grad_phi_variable = var_db->getVariable("grad phi (LSMLIB)");
  } else {
    grad_phi_variable = new CellVariable<DIM,LSMLIB_REAL>("grad phi (LSMLIB)",DIM); 
  }
 
  // reserve memory for scratch variable PatchData Handles
  d_phi_handles.reserve(d_tvd_runge_kutta_order);

  // phi - "CURRENT" context for time advance
  d_phi_handles[0] = var_db->registerVariableAndContext(
      phi_variable, current_context, d_level_set_ghostcell_width);
  d_solution_variables.setFlag(d_phi_handles[0]);

  // phi - "SCRATCH" context for time advance
  for (int k=1; k < d_tvd_runge_kutta_order; k++) {
    stringstream context_name("");
    context_name << "TVD_RK_SCRATCH_" << k;
    d_phi_handles[k] = var_db->registerVariableAndContext(
      phi_variable, 
      var_db->getContext(context_name.str()),
      d_level_set_ghostcell_width);
    d_time_advance_scratch_variables.setFlag(d_phi_handles[k]);
  }

  // upwind grad(phi)
  d_grad_phi_upwind_handle = var_db->registerVariableAndContext(
    grad_phi_variable, upwind_context, zero_ghostcell_width);
  if (d_lsm_velocity_field_strategy->providesExternalVelocityField()) {
    d_time_advance_scratch_variables.setFlag(d_grad_phi_upwind_handle);
  }

  // forward and backward grad(phi)
  d_grad_phi_plus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, plus_context, zero_ghostcell_width);
  d_grad_phi_minus_handle = var_db->registerVariableAndContext(
    grad_phi_variable, minus_context, zero_ghostcell_width);
  if (d_lsm_velocity_field_strategy->providesNormalVelocityField()) {
    d_compute_stable_dt_scratch_variables.setFlag(d_grad_phi_plus_handle);
    d_compute_stable_dt_scratch_variables.setFlag(d_grad_phi_minus_handle);
    d_time_advance_scratch_variables.setFlag(d_grad_phi_plus_handle);
    d_time_advance_scratch_variables.setFlag(d_grad_phi_minus_handle);
  } else {
    d_reinitialization_scratch_variables.setFlag(d_grad_phi_plus_handle);
    d_reinitialization_scratch_variables.setFlag(d_grad_phi_minus_handle);
    d_orthogonalization_scratch_variables.setFlag(d_grad_phi_plus_handle);
    d_orthogonalization_scratch_variables.setFlag(d_grad_phi_minus_handle);
  }

  // RHS for phi updates
  Pointer< CellVariable<DIM,LSMLIB_REAL> > rhs_phi_variable;
  if (var_db->checkVariableExists("rhs phi (LSMLIB)")) {
    rhs_phi_variable = var_db->getVariable("rhs phi (LSMLIB)");
  } else {
    rhs_phi_variable = new CellVariable<DIM,LSMLIB_REAL>("rhs phi (LSMLIB)",1); 
  }
  d_rhs_phi_handle = var_db->registerVariableAndContext(
    rhs_phi_variable, scratch_context, zero_ghostcell_width);
  d_time_advance_scratch_variables.setFlag(d_rhs_phi_handle);


  /* 
   * Initialize psi variables for codimension-two problems
   */
  // reserve memory for scratch variable PatchData Handles
  d_psi_handles.reserve(d_tvd_runge_kutta_order);

  if (d_codimension == 2) {

    // create variables
    Pointer< CellVariable<DIM,LSMLIB_REAL> > psi_variable;
    if (var_db->checkVariableExists("psi (LSMLIB)")) {
      psi_variable = var_db->getVariable("psi (LSMLIB)");
    } else {
      psi_variable = new CellVariable<DIM,LSMLIB_REAL>(
        "psi (LSMLIB)", d_num_level_set_fcn_components); 
    }
    Pointer< CellVariable<DIM,LSMLIB_REAL> > grad_psi_variable;
    if (var_db->checkVariableExists("grad psi (LSMLIB)")) {
      grad_psi_variable = var_db->getVariable("grad psi (LSMLIB)");
    } else {
      grad_psi_variable = new CellVariable<DIM,LSMLIB_REAL>(
        "grad psi (LSMLIB)",DIM); 
    } 

    // psi - "CURRENT" context for time advance
    d_psi_handles[0] = var_db->registerVariableAndContext(
        psi_variable, current_context, d_level_set_ghostcell_width);
    d_solution_variables.setFlag(d_psi_handles[0]);

    // psi - "SCRATCH" context for time advance
    for (int k=1; k < d_tvd_runge_kutta_order; k++) {
      stringstream context_name("");
      context_name << "TVD_RK_SCRATCH_" << k;
      d_psi_handles[k] = var_db->registerVariableAndContext(
        psi_variable, 
        var_db->getContext(context_name.str()),
        d_level_set_ghostcell_width);
      d_time_advance_scratch_variables.setFlag(d_psi_handles[k]); 
    }

    // upwind grad(psi)
    d_grad_psi_upwind_handle = var_db->registerVariableAndContext(
      grad_psi_variable, upwind_context, zero_ghostcell_width);
    if (d_lsm_velocity_field_strategy->providesExternalVelocityField()) {
      d_time_advance_scratch_variables.setFlag(d_grad_psi_upwind_handle);
    }

    // forward and backward grad(phi)
    d_grad_psi_plus_handle = var_db->registerVariableAndContext(
      grad_psi_variable, plus_context, zero_ghostcell_width);
    d_grad_psi_minus_handle = var_db->registerVariableAndContext(
      grad_psi_variable, minus_context, zero_ghostcell_width);
    if (d_lsm_velocity_field_strategy->providesNormalVelocityField()) {
      d_compute_stable_dt_scratch_variables.setFlag(d_grad_psi_plus_handle);
      d_compute_stable_dt_scratch_variables.setFlag(d_grad_psi_minus_handle);
      d_time_advance_scratch_variables.setFlag(d_grad_psi_plus_handle);
      d_time_advance_scratch_variables.setFlag(d_grad_psi_minus_handle);
    } else {
      d_reinitialization_scratch_variables.setFlag(d_grad_psi_plus_handle);
      d_reinitialization_scratch_variables.setFlag(d_grad_psi_minus_handle);
      d_orthogonalization_scratch_variables.setFlag(d_grad_psi_plus_handle);
      d_orthogonalization_scratch_variables.setFlag(d_grad_psi_minus_handle);
    }

    // RHS for psi updates
    Pointer< CellVariable<DIM,LSMLIB_REAL> > rhs_psi_variable;
    if (var_db->checkVariableExists("rhs psi (LSMLIB)")) {
      rhs_psi_variable = var_db->getVariable("rhs psi (LSMLIB)");
    } else {
      rhs_psi_variable = new CellVariable<DIM,LSMLIB_REAL>("rhs psi (LSMLIB)",1); 
    }
    d_rhs_psi_handle = var_db->registerVariableAndContext(
      rhs_psi_variable, scratch_context, zero_ghostcell_width);
    d_time_advance_scratch_variables.setFlag(d_rhs_psi_handle);

  } else { // set PatchData handles for filling psi scratch data to -1 
           // (a bogus value)
    for (int k=0; k < d_tvd_runge_kutta_order; k++) 
      d_psi_handles[k] = -1;
  }

  /*
   * Initialize control volume variables
   */
  Pointer< CellVariable<DIM,LSMLIB_REAL> > control_volume;
  if (var_db->checkVariableExists("control volume (LSMLIB)")) {
    control_volume = var_db->getVariable("control volume (LSMLIB)");
  } else {
    control_volume = new CellVariable<DIM,LSMLIB_REAL>("control volume (LSMLIB)",1);
  }
  d_control_volume_handle = var_db->registerVariableAndContext(
    control_volume, current_context, zero_ghostcell_width);
  d_persistent_variables.setFlag(d_control_volume_handle);

  /*
   * Register phi, psi, and control volume as restart PatchData items.
   */
  var_db->registerPatchDataForRestart(d_phi_handles[0]);
  if (d_codimension == 2) {
    var_db->registerPatchDataForRestart(d_psi_handles[0]);
  }
  var_db->registerPatchDataForRestart(d_control_volume_handle);

}

/* initializeCommunicationObjects() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::initializeCommunicationObjects()
{
  // lookup refine operations
  Pointer< RefineOperator<DIM> > refine_op =
    d_grid_geometry->lookupRefineOperator(
      VariableDatabase<DIM>::getDatabase()->getVariable("phi (LSMLIB)"),
      "LINEAR_REFINE");

  // set up communications objects for filling a new level (used during
  // initialization of a level and in initializing velocity fields)
  d_fill_new_level = new RefineAlgorithm<DIM>;

  // fill algorithm for initializing new level
  d_fill_new_level->registerRefine(
    d_phi_handles[0], d_phi_handles[0], d_phi_handles[0], 
    refine_op);
  if (d_codimension == 2) {
    d_fill_new_level->registerRefine(
      d_psi_handles[0], d_psi_handles[0], d_psi_handles[0], 
      refine_op);
  }

  // set up objects for filling boundary data during the calculation
  // of the stable time step 
  d_fill_bdry_compute_stable_dt = new RefineAlgorithm<DIM>;

  // empty out the boundary bdry fill schedules 
  d_fill_bdry_sched_compute_stable_dt.setNull();

  d_fill_bdry_compute_stable_dt->registerRefine(
    d_phi_handles[0], d_phi_handles[0], 
    d_phi_handles[0], refine_op);
  if (d_codimension == 2) {
    d_fill_bdry_compute_stable_dt->registerRefine(
      d_psi_handles[0], d_psi_handles[0], 
      d_psi_handles[0], refine_op);
  }

  // set up objects for filling boundary data during the 
  // time advance of the level set functions
  d_fill_bdry_time_advance.resizeArray(d_tvd_runge_kutta_order);
  d_fill_bdry_sched_time_advance.resizeArray(d_tvd_runge_kutta_order);

  for (int k = 0; k < d_tvd_runge_kutta_order; k++) {
    d_fill_bdry_time_advance[k] = new RefineAlgorithm<DIM>;

    // empty out the boundary bdry fill schedules 
    d_fill_bdry_sched_time_advance[k].setNull();

    // fill algorithm before time advance
    d_fill_bdry_time_advance[k]->registerRefine(
      d_phi_handles[k], 
      d_phi_handles[k], 
      d_phi_handles[k], 
      refine_op);
    if (d_codimension == 2) {
      d_fill_bdry_time_advance[k]->registerRefine(
        d_psi_handles[k], 
        d_psi_handles[k], 
        d_psi_handles[k], 
        refine_op);
    }  

  } // end loop setting up data transfers for TVD Runge-Kutta time advance

}


/* getFromInput() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::getFromInput(
  Pointer<Database> db, 
  bool is_from_restart)
{
  /*
   * Read in parameters that override values in the restart file
   */
  if (is_from_restart) {
    if (db->keyExists("end_time")) d_end_time = db->getDouble("end_time");
  } else {
    d_end_time = db->getDoubleWithDefault("end_time", LSM_DEFAULT_END_TIME);
  }

  if (is_from_restart) {
    if (db->keyExists("reinitialization_interval")) {
      d_reinitialization_interval = 
        db->getInteger("reinitialization_interval");
      d_use_reinitialization = (d_reinitialization_interval > 0);
    }
    if (db->keyExists("orthogonalization_interval")) {
      d_orthogonalization_interval = 
        db->getInteger("orthogonalization_interval");
      d_use_orthogonalization = ( (d_orthogonalization_interval > 0) &&
                                  (d_codimension == 2) );
    }
  } else {
    d_reinitialization_interval = db->getIntegerWithDefault(
      "reinitialization_interval", LSM_DEFAULT_REINITIALIZATION_INTERVAL);
    d_use_reinitialization = (d_reinitialization_interval > 0);
    d_orthogonalization_interval = db->getIntegerWithDefault(
      "orthogonalization_interval", LSM_DEFAULT_ORTHOGONALIZATION_INTERVAL);
    d_use_orthogonalization = ( (d_orthogonalization_interval > 0) &&
                                (d_codimension == 2) );
  }

  if (is_from_restart) {

    if (db->keyExists("reinitialization_stop_tol")) {
      d_reinitialization_stop_tol = db->getDouble("reinitialization_stop_tol");
      d_use_reinitialization_stop_tol = true;
    } 
    if (db->keyExists("reinitialization_stop_dist")) {
      d_reinitialization_stop_dist = 
        db->getDouble("reinitialization_stop_dist");
      d_use_reinitialization_stop_dist = true;
    } 
    if (db->keyExists("reinitialization_max_iters")) {
      d_reinitialization_max_iters = 
        db->getInteger("reinitialization_max_iters");
      d_use_reinitialization_max_iters = true;
    } 
    if (db->keyExists("orthogonalization_stop_tol")) {
      d_orthogonalization_stop_tol = 
        db->getDouble("orthogonalization_stop_tol");
      d_use_orthogonalization_stop_tol = true;
    } 
    if (db->keyExists("orthogonalization_stop_dist")) {
      d_orthogonalization_stop_dist = 
        db->getDouble("orthogonalization_stop_dist");
      d_use_orthogonalization_stop_dist = true;
    } 
    if (db->keyExists("orthogonalization_max_iters")) {
      d_orthogonalization_max_iters = 
        db->getInteger("orthogonalization_max_iters");
      d_use_orthogonalization_max_iters = true;
    } 

  } else {

    if (db->keyExists("reinitialization_stop_tol")) {
      d_reinitialization_stop_tol = db->getDouble("reinitialization_stop_tol");
      d_use_reinitialization_stop_tol = true;
    } else {
      d_reinitialization_stop_tol = 0.0;
      d_use_reinitialization_stop_tol = false;
    }
    if (db->keyExists("reinitialization_stop_dist")) {
      d_reinitialization_stop_dist = 
       db->getDouble("reinitialization_stop_dist");
      d_use_reinitialization_stop_dist = true;
    } else {
      d_reinitialization_stop_dist = 0.0;
      d_use_reinitialization_stop_dist = false;
    }
    if (db->keyExists("reinitialization_max_iters")) {
      d_reinitialization_max_iters = 
       db->getInteger("reinitialization_max_iters");
      d_use_reinitialization_max_iters = true;
    } else {
      d_reinitialization_max_iters = 0;
      d_use_reinitialization_max_iters = false;
    }
    if (db->keyExists("orthogonalization_stop_tol")) {
      d_orthogonalization_stop_tol = 
       db->getDouble("orthogonalization_stop_tol");
      d_use_orthogonalization_stop_tol = true;
    } else {
      d_orthogonalization_stop_tol = 0.0;
      d_use_orthogonalization_stop_tol = false;
    }
    if (db->keyExists("orthogonalization_stop_dist")) {
      d_orthogonalization_stop_dist = 
        db->getDouble("orthogonalization_stop_dist");
      d_use_orthogonalization_stop_dist = true;
    } else {
      d_orthogonalization_stop_dist = 0.0;
      d_use_orthogonalization_stop_dist = false;
    }
    if (db->keyExists("orthogonalization_max_iters")) {
      d_orthogonalization_max_iters = 
        db->getInteger("orthogonalization_max_iters");
      d_use_orthogonalization_max_iters = true;
    } else {
      d_orthogonalization_max_iters = 0;
      d_use_orthogonalization_max_iters = false;
    }

  }  // not from restart case case for reinit/ortho parameters
  
  // if no stopping criteria were specified use the default number of
  // maximum iterations
  if ( !( d_use_reinitialization_stop_tol || 
          d_use_reinitialization_stop_dist ||
          d_use_reinitialization_max_iters) ) {
    d_use_reinitialization_max_iters = true;
    d_reinitialization_max_iters = LSM_DEFAULT_REINITIALIZATION_MAX_ITERS;
  }
  if ( !( d_use_orthogonalization_stop_tol || 
          d_use_orthogonalization_stop_dist ||
          d_use_orthogonalization_max_iters) ) {
    d_use_orthogonalization_max_iters = true;
    d_orthogonalization_max_iters = LSM_DEFAULT_ORTHOGONALIZATION_MAX_ITERS;
  }

  // get verbose mode
  if (is_from_restart) {
    if (db->keyExists("verbose_mode")) d_verbose_mode = 
      db->getBool("verbose_mode");
  } else {
    d_verbose_mode = db->getBoolWithDefault("verbose_mode", 
      LSM_DEFAULT_VERBOSE_MODE);
  } 


  /*
   * If computation is NOT from restart, read in all of the 
   * user-specified parameters from the input database.
   */
  if (!is_from_restart) {
    // read in level set parameters
    d_start_time = db->getDoubleWithDefault("start_time", 
                                            LSM_DEFAULT_START_TIME);
    d_cfl_number = db->getDoubleWithDefault("cfl_number", 
                                            LSM_DEFAULT_CFL_NUMBER);
  
    string spatial_derivative_type = db->getStringWithDefault(
      "spatial_derivative_type", LSM_DEFAULT_SPATIAL_DERIVATIVE_TYPE);
    if (spatial_derivative_type == "WENO") {
      d_spatial_derivative_type = WENO;
      d_spatial_derivative_order = db->getIntegerWithDefault(
        "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_WENO_ORDER);
    } else if (spatial_derivative_type == "ENO") {
      d_spatial_derivative_type = ENO;
      d_spatial_derivative_order = db->getIntegerWithDefault(
        "spatial_derivative_order", LSM_DEFAULT_SPATIAL_DERIVATIVE_ENO_ORDER);
    } else {
      d_spatial_derivative_type = UNKNOWN;
    }
    d_tvd_runge_kutta_order = db->getIntegerWithDefault(
      "tvd_runge_kutta_order", LSM_DEFAULT_TVD_RUNGE_KUTTA_ORDER);

    // check that spatial derivative type, spatial derivative order,
    // and TVD Runge-Kutta order are valid.
    if (d_spatial_derivative_type == ENO) {
      if ( (d_spatial_derivative_order < 1) ||
           (d_spatial_derivative_order > 3) ) {
        TBOX_ERROR(d_object_name 
                << "::getFromInput(): "
                << "Unsupported order for ENO derivative.  "
                << "Only ENO1, ENO2, and ENO3 supported."
                << endl );
      }
    } else if (d_spatial_derivative_type == WENO) {
      if ( (d_spatial_derivative_order != 5) )
        TBOX_ERROR(d_object_name 
                << "::getFromInput(): "
                << "Unsupported order for WENO derivative.  "
                << "Only WENO5 supported."
                << endl );
    }
    if ( (d_tvd_runge_kutta_order < 1) ||
         (d_tvd_runge_kutta_order > 3) ) {
      TBOX_ERROR(d_object_name 
              << "::getFromInput(): "
              << "Unsupported TVD Runge-Kutta order.  "
              << "Only TVD-RK1, TVD-RK2, and TVD-RK3 supported."
              << endl );
    }
  
    // read in boundary conditions
    for (int comp=0; comp < d_num_level_set_fcn_components; comp++) {
      stringstream lower_bc_phi_key;
      lower_bc_phi_key << "lower_bc_phi_";
      lower_bc_phi_key << comp; 
      int lower_bc[DIM];
      if (db->keyExists(lower_bc_phi_key.str())) {
        db->getIntegerArray(lower_bc_phi_key.str(), lower_bc, DIM);
        for (int dim = 0; dim < DIM; dim++) {
          d_lower_bc_phi[comp](dim) = lower_bc[dim];
        }
      }

      stringstream upper_bc_phi_key;
      upper_bc_phi_key << "upper_bc_phi_";
      upper_bc_phi_key << comp; 
      int upper_bc[DIM];
      if (db->keyExists(upper_bc_phi_key.str())) {
        db->getIntegerArray(upper_bc_phi_key.str(), upper_bc, DIM);
        for (int dim = 0; dim < DIM; dim++) {
          d_upper_bc_phi[comp](dim) = upper_bc[dim];
        }
      }

      stringstream lower_bc_psi_key;
      lower_bc_psi_key << "lower_bc_psi_";
      lower_bc_psi_key << comp; 
      if (db->keyExists(lower_bc_psi_key.str())) {
        db->getIntegerArray(lower_bc_psi_key.str(), lower_bc, DIM);
        for (int dim = 0; dim < DIM; dim++) {
          d_lower_bc_psi[comp](dim) = lower_bc[dim];
        }
      }

      stringstream upper_bc_psi_key;
      upper_bc_psi_key << "upper_bc_psi_";
      upper_bc_psi_key << comp; 
      if (db->keyExists(upper_bc_psi_key.str())) {
        db->getIntegerArray(upper_bc_psi_key.str(), upper_bc, DIM);
        for (int dim = 0; dim < DIM; dim++) {
          d_upper_bc_psi[comp](dim) = upper_bc[dim];
        }
      }
    } // end: read boundary conditions


    // read in AMR parameters
    d_use_AMR = db->getBoolWithDefault("use_AMR", LSM_DEFAULT_USE_AMR);
    d_regrid_interval = db->getIntegerWithDefault("regrid_interval", 
      LSM_DEFAULT_REGRID_INTERVAL);
    d_tag_buffer_width = db->getIntegerWithDefault("tag_buffer_width", 
      LSM_DEFAULT_TAG_BUFFER_WIDTH);
    d_refinement_cutoff_value = db->getDoubleWithDefault(
      "refinement_cutoff_value", LSM_DEFAULT_REFINEMENT_CUTOFF_VALUE);
  } // end case (NOT FROM RESTART)

}


/* getFromRestart() */
template <int DIM> 
void LevelSetFunctionIntegrator<DIM>::getFromRestart()
{

  // open restart file
  Pointer<Database> root_db =
    RestartManager::getManager()->getRootDatabase();

  Pointer<Database> db;
  if ( root_db->isDatabase(d_object_name) ) {
    db = root_db->getDatabase(d_object_name);
  } else {
    TBOX_ERROR(d_object_name 
            << "::getFromRestart(): "
            << "Restart database corresponding to "
            << d_object_name << " not found in restart file" 
            << endl);
  }

  // check that version of restart file is consistent with 
  // version of class
  int ver = db->getInteger("LEVEL_SET_METHOD_INTEGRATOR_VERSION");
  if (ver != LEVEL_SET_METHOD_INTEGRATOR_VERSION) {
    TBOX_ERROR(d_object_name 
            << "::getFromRestart(): "
            << "Restart file version different than class version." 
            << endl);
  }

  /*
   * Read in user-specified parameters 
   */

  // read in level set parameters
  d_start_time = db->getDouble("d_start_time");
  d_end_time = db->getDouble("d_end_time");
  d_cfl_number = db->getDouble("d_cfl_number");

  d_spatial_derivative_type = 
    (SPATIAL_DERIVATIVE_TYPE) db->getInteger("d_spatial_derivative_type");
  d_spatial_derivative_order = db->getInteger("d_spatial_derivative_order");
  d_tvd_runge_kutta_order = db->getInteger("d_tvd_runge_kutta_order");

  d_reinitialization_interval = db->getInteger("d_reinitialization_interval");
  d_reinitialization_stop_tol = db->getDouble("d_reinitialization_stop_tol");
  d_reinitialization_stop_dist = db->getDouble("d_reinitialization_stop_dist");
  d_reinitialization_max_iters= db->getInteger("d_reinitialization_max_iters");
  d_orthogonalization_interval = db->getInteger("d_orthogonalization_interval");
  d_orthogonalization_stop_tol = db->getDouble("d_orthogonalization_stop_tol");
  d_orthogonalization_stop_dist = 
    db->getDouble("d_orthogonalization_stop_dist");
  d_orthogonalization_max_iters = 
    db->getInteger("d_orthogonalization_max_iters");

  d_use_AMR = db->getBool("d_use_AMR");
  d_regrid_interval = db->getInteger("d_regrid_interval");
  d_tag_buffer_width = db->getInteger("d_tag_buffer_width");
  d_refinement_cutoff_value = db->getDouble("d_refinement_cutoff_value");

  d_verbose_mode = db->getBool("d_verbose_mode");

  /*
   * Read in state parameters
   */
  d_current_time = db->getDouble("d_current_time");
  d_num_integration_steps_taken = db->getInteger("d_num_integration_steps_taken");
  d_reinitialization_count = db->getInteger("d_reinitialization_count");
  d_orthogonalization_count = db->getInteger("d_orthogonalization_count");
  d_regrid_count = db->getInteger("d_regrid_count");

  d_use_reinitialization = db->getBool("d_use_reinitialization");
  d_use_reinitialization_stop_tol = 
    db->getBool("d_use_reinitialization_stop_tol");
  d_use_reinitialization_stop_dist =
    db->getBool("d_use_reinitialization_stop_dist");
  d_use_reinitialization_max_iters = 
    db->getBool("d_use_reinitialization_max_iters");
  d_use_orthogonalization = db->getBool("d_use_orthogonalization");
  d_use_orthogonalization_stop_tol = 
    db->getBool("d_use_orthogonalization_stop_tol");
  d_use_orthogonalization_stop_dist =
    db->getBool("d_use_orthogonalization_stop_dist");
  d_use_orthogonalization_max_iters =
    db->getBool("d_use_orthogonalization_max_iters");
}

} // end LSMLIB namespace

#endif
