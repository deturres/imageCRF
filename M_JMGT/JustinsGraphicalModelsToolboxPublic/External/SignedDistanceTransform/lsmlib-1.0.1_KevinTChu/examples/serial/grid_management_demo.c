/*
 * File:        grid_management_demo.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision: 149 $
 * Modified:    $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
 * Description: Demo code for grid management functions.
 */

#include <stdio.h>
#include "LSMLIB_config.h"
#include "lsm_file.h"
#include "lsm_grid.h"

int main(void)
{
  Grid *g_original, *g_from_ascii_file, *g_from_binary_file;
  int num_dims;
  int grid_dims[3];
  LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy_type;
  LSMLIB_REAL x_lo[3];
  LSMLIB_REAL x_hi[3];


  printf("\n************ Demo 2D Grid management functions ************\n\n");

  num_dims = 2;
  grid_dims[0] = 10; grid_dims[1] = 20;
  accuracy_type = VERY_HIGH;
  x_lo[0] = -1; x_lo[1] = -1;
  x_hi[0] = 1; x_hi[1] = 1;

  g_original = createGridSetGridDims(num_dims, grid_dims, 
                                     x_lo, x_hi, accuracy_type);
  printf("====================== Original Grid =========================\n");
  printGrid(g_original,stdout); 
  printf("==============================================================\n");

  writeGridToAsciiFile(g_original,"grid_2d_demo.ascii",NO_ZIP);
  g_from_ascii_file = readGridFromAsciiFile("grid_2d_demo.ascii");
  printf("=================== Grid From ASCII File =====================\n");
  printGrid(g_from_ascii_file,stdout);
  printf("==============================================================\n");

  writeGridToBinaryFile(g_original,"grid_2d_demo.binary",GZIP);
  g_from_binary_file = readGridFromBinaryFile("grid_2d_demo.binary");
  printf("=================== Grid From Binary File ====================\n");
  printGrid(g_from_binary_file,stdout);
  printf("==============================================================\n");

  printf("\n**************** End 2D Grid management demo ****************\n\n");


  printf("\n************ Demo 3D Grid management functions ************\n\n");

  num_dims = 3;
  grid_dims[0] = 10; grid_dims[1] = 20; grid_dims[2] = 30;
  accuracy_type = VERY_HIGH;
  x_lo[0] = -1; x_lo[1] = -1; x_lo[2] = -1;
  x_hi[0] = 1; x_hi[1] = 1; x_hi[2] = 1;

  g_original = createGridSetGridDims(num_dims, grid_dims, 
                                     x_lo, x_hi, accuracy_type);
  printf("====================== Original Grid =========================\n");
  printGrid(g_original,stdout); 
  printf("==============================================================\n");

  writeGridToAsciiFile(g_original,"grid_3d_demo.ascii",NO_ZIP);
  g_from_ascii_file = readGridFromAsciiFile("grid_3d_demo.ascii");
  printf("=================== Grid From ASCII File =====================\n");
  printGrid(g_from_ascii_file,stdout);
  printf("==============================================================\n");

  writeGridToBinaryFile(g_original,"grid_3d_demo.binary",GZIP);
  g_from_binary_file = readGridFromBinaryFile("grid_3d_demo.binary");
  printf("=================== Grid From Binary File ====================\n");
  printGrid(g_from_binary_file,stdout);
  printf("==============================================================\n");

  printf("\n**************** End 3D Grid management demo ****************\n\n");

  return 0;
}
