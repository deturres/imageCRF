/*
 * File:        Array__Pointer__RefineAlgorithm.NDIM.cc
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Explicit template instantiation of LSMLIB classes 
 */

#include "SAMRAI_config.h"
#include "RefineAlgorithm.h"
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"

namespace SAMRAI {

template class tbox::Array< tbox::Pointer< 
  xfer::RefineAlgorithm<NDIM> 
> >;

}
