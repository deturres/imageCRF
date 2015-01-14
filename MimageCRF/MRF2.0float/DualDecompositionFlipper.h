/*
 * DualDecompositionFlipper.h
 *
 *  Created on: May 6, 2012
 *      Author: bhole
 */

#ifndef DUALDECOMPOSITIONFLIPPER_H_
#define DUALDECOMPOSITIONFLIPPER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "DualDecomposition.h"

#define FloatType float
#define FLOATTYPE float

enum SplitType { SUBNONSUB, TREESADD };

class DualDecompositionFlipper : public DualDecomposition {
public:
  typedef CostVal REAL;
  DualDecompositionFlipper(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecompositionFlipper(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecompositionFlipper();

  void setDDFlipperFunction(MRF::SmoothCostGeneralFn fnCost_dd_flipper);
  void setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper);
  void setNodeMappingStructure(std::vector<int> **sg_node_to_g_node);
  void setflippedinfo(std::vector<int> **flippedinfo);


  // subnonsub breaks graph into one submodular graph and minimum number of nonsub subgraphs
  // for a planar grid graph, only one nonsubmod graph is sufficient
  // TREESADD starts with the tree splitting used by DD-T and then converts them to flipper graphs

  void setSplitType(SplitType s_var);

protected:

  virtual void optimizeAlg(int nIterations);

  std::vector<int> **subgraph_node_to_graph_node_ptop;
  SmoothCostGeneralFn m_ddflippersmoothFn;
  SmoothCostGeneralFn m_flippersmoothFn;
  std::vector<int> **flippedinfo_ptop;

private:

  SplitType m_splittype;

  void createSubgraphs(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<int> > *subgraphs_edge_list,
      std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes);

  void createSubgraphsTreeAdd(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<int> > *subgraphs_edge_list,
      std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes);


  void moveTreeToSubgraph(
      std::vector<int> *tree_edge_list,
      std::vector<std::vector<int> > *subgraphs_edge_list);

  void markFlippedNodesFlipper(
      const std::vector<int> &subgraph_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &edge_type,
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      std::vector<int> *mark_flip_nodes);

  void traverseGraph(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<int> &subgraphs_edge_list0,
      const std::vector<std::vector<int> > &edge_list);

  void moveTreeEdgesToFlipper(
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &edge_type,
      std::vector<int> *mark_flip_nodes,
      std::vector<int> *subgraph_edge_list,
      std::vector<int> *tree_edge_list);

  void computeDataCostForSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
      std::vector<MRF::CostVal*> *D_a,
      std::vector<std::vector<int> > *subgraph_node_to_graph_node,
      std::vector<DataCost*> *dcost_vec);

  void getSubGraphNodeNo(const std::vector<std::vector<int> > &node_in_subgraphs_g,
      const int &subgraph_no, int *subgraph_node1);

  void dumpSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
                 const std::vector<std::vector<int> > &edge_list);

  void displayNodeEdgeType(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &edge_type);

  void getWhatSubGraphsNodeIsUsedIn(
      const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
      std::vector<int> *num_nodes_in_subgraphs);


  /*
  void extractSubGraphEdges(const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *global_used_edges,
      std::vector<std::vector<int> > *subgraphs_edge_list,
      std::vector<int> *subgraph_type);

  void markFlippedNodes(const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &subgraph_type,
      std::vector<std::vector<std::vector<int> > >*node_in_subgraphs);


*/

};



#endif /* DUALDECOMPOSITIONFLIPPER_H_ */
