/*
 * File:        PatchModule.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Implementation for concrete subclass of 
 *              LevelSetMethodPatchStrategy that computes the single patch
 *              numerical routines for the 2d level set method example problem
 */


#include "PatchModule.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"

// headers for level set method numerical kernels
extern "C" {
  #include "patchmodule_fort.h"
}

// SAMRAI namespaces
using namespace geom; 
using namespace pdat; 
using namespace LSMLIB;

// CONSTANTS
const LSMLIB_REAL PatchModule::s_default_radius = 0.25;

PatchModule::PatchModule(
  Pointer<Database> input_db,
  const string& object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!input_db.isNull());
  assert(!object_name.empty());
#endif

  // set object name and grid geometry
  d_object_name = object_name;

  // read in input data
  getFromInput(input_db);

}

void PatchModule::initializeLevelSetFunctionsOnPatch(
  Patch<2>& patch,
  const LSMLIB_REAL data_time,
  const int phi_handle,
  const int psi_handle)
{
  Pointer< CellData<2,LSMLIB_REAL> > level_set_data =
    patch.getPatchData( phi_handle );

  LSMLIB_REAL* level_set_data_ptr = level_set_data->getPointer();

  Pointer< CartesianPatchGeometry<2> > patch_geom 
    = patch.getPatchGeometry();

#ifdef LSMLIB_DOUBLE_PRECISION
  const double* dx = patch_geom->getDx();
  const double* x_lower = patch_geom->getXLower();
#else
  const double* dx_double = patch_geom->getDx();
  const double* x_lower_double = patch_geom->getXLower();
  float dx[2], x_lower[2];
  dx[0] = dx_double[0]; dx[1] = dx_double[1];
  x_lower[0] = x_lower_double[0]; x_lower[1] = x_lower_double[1];
#endif

  Box<2> box = level_set_data->getBox();
  Box<2> ghostbox = level_set_data->getGhostBox();
  const IntVector<2> box_lower = box.lower();
  const IntVector<2> box_upper = box.upper();
  const IntVector<2> ghostbox_lower = ghostbox.lower();
  const IntVector<2> ghostbox_upper = ghostbox.upper();

  switch (d_initial_level_set) 
  {
    case CIRCLE: // circle
    {
      INIT_CIRCLE(
        level_set_data_ptr,
        &ghostbox_lower[0],
        &ghostbox_upper[0],
        &ghostbox_lower[1],
        &ghostbox_upper[1],
        &box_lower[0],
        &box_upper[0],
        &box_lower[1],
        &box_upper[1],
        x_lower,
        dx,
        d_center,
        &d_radius);
      break;
    }
    default: {}
  }; 

}

void PatchModule::setLevelSetFunctionBoundaryConditions(
    Patch<2>& patch,
    const LSMLIB_REAL fill_time,
    const int phi_handle,
    const int psi_handle,
    const IntVector<2>& ghost_width_to_fill)
{
}

void PatchModule::printClassData(ostream &os) const
{
  os << "\nPatchModule::printClassData..." << endl;
  os << "PatchModule: this = " << (PatchModule*)this 
     << endl;
  os << "d_object_name = " << d_object_name << endl;
  os << "d_initial_level_set = " << d_initial_level_set << endl;

  // KTC - PUT MORE HERE
  os << endl;
}


void PatchModule::getFromInput(
  Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif

  // set initial level_set_selector
  d_initial_level_set = db->getIntegerWithDefault("initial_level_set", 0);

  // get auxilliary parameters for initial level set
  switch (d_initial_level_set) {
    case CIRCLE: {

#ifdef LSMLIB_DOUBLE_PRECISION
      d_radius = db->getDoubleWithDefault("radius", s_default_radius);
#else
      d_radius = db->getFloatWithDefault("radius", s_default_radius);
#endif

      if (db->keyExists("center")) {
#ifdef LSMLIB_DOUBLE_PRECISION
        db->getDoubleArray("center", d_center, 2);
#else
        db->getFloatArray("center", d_center, 2);
#endif

      } else {
        d_center[0] = 0.0;
        d_center[1] = 0.0;
      }
      break;
    }

    default: {
      TBOX_ERROR(  "PatchModule"
                << "::getFromInput()"
                << ":Invalid type of initial level set"
                << endl );
    }
  };

}

