/*
 * FlipperForDD.h
 *
 *  Created on: May 18, 2012
 *      Author: bhole
 */

#ifndef FLIPPERFORDD_H_
#define FLIPPERFORDD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "MRFWrapper.h"


#define FloatType float
#define FLOATTYPE float


class FlipperForDD : public MRFWrapper {
public:
  typedef CostVal REAL;
  FlipperForDD(int width, int height, int nLabels, EnergyFunction *eng);
  FlipperForDD(int nPixels, int nLabels,EnergyFunction *eng);
  ~FlipperForDD();

  void setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper);
  void setflippedinfo(std::vector<int> **flippedinfo_f);

  void setUpFlipperForDD();

  void updateFlipperDataStructures(const int idx, const int pixelno, const int d,
      const float value);

protected:

  virtual void optimizeAlg(int nIterations);
  SmoothCostGeneralFn m_flippersmoothFn;
  std::vector<int> **flippedinfo_ptop;

private:

  void extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
                       std::vector<std::vector<int> > *edge_list);

  int checkValidFlipperGraph();


  void markFlippedNodes(
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


  // [node_no, [other_no, edge_no]*]*
  std::vector<std::vector<std::vector<int> > > node_edge_list;

  // [edge_no, [first_node, second_node]]*
  std::vector<std::vector<int> > edge_list;


  // 0 indicates not used, 1 indicates edge is used by at least one tree.
  //std::vector<int> global_used_edges(edge_list.size(), 0);

  // ie [edge_no [type]]* - all edges included in flipper graph
  // type is 0 if submodular and 1 for nonsubmodular
  std::vector<int>  edge_list_type;

  // stores the nodes that are flipped as 1 else 0
  std::vector<int> nodes_flipped;

  DataCost* dcost_flipper;
  EnergyFunction* energy_flipper;
  MRF* mrf_flipper;
  MRF::CostVal* D_a;

  SmoothnessCost* sscost;

};




#endif /* FLIPPERFORDD_H_ */
