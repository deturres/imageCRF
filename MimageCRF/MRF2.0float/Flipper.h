/*
 * Flipper.h
 *
 *  Created on: May 6, 2012
 *      Author: bhole
 */

#ifndef FLIPPER_H_
#define FLIPPER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "MRFWrapper.h"


#define FloatType float
#define FLOATTYPE float


class Flipper : public MRFWrapper {
public:
  typedef CostVal REAL;
  Flipper(int width, int height, int nLabels, EnergyFunction *eng);
  Flipper(int nPixels, int nLabels,EnergyFunction *eng);
  ~Flipper();

  void setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper);
  void setflippedinfo(std::vector<int> **flippedinfo_f);

  // void updatenodeArray(int pixelno, int d, float value);

protected:

  virtual void optimizeAlg(int nIterations);
  SmoothCostGeneralFn m_flippersmoothFn;
  std::vector<int> **flippedinfo_ptop;

private:

  void extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
                       std::vector<std::vector<int> > *edge_list);

  int checkValidFlipperGraph();


  int markFlippedNodes(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<int> &edge_list_type,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<int> *nodes_flipped);

  void computeDataCostForFlipper(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<int> &edge_list_type,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &nodes_flipped,
      MRF::CostVal* D_a);

};


#endif /* FLIPPER_H_ */
