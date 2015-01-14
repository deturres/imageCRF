/*
 * DualDecompositionFlipperES.h
 *
 *  Created on: May 29, 2012
 *      Author: bhole
 */

// edge sharing enabled.


#ifndef DUALDECOMPOSITIONFLIPPERES_H_
#define DUALDECOMPOSITIONFLIPPERES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "common_typedefs.h"
#include "DualDecomposition.h"

#define FloatType float
#define FLOATTYPE float

namespace DDFES {
enum SplitType { SUBNONSUB, TREESADD, TREESADD_EDGESHARE };
}

class DualDecompositionFlipperES : public DualDecomposition {
public:
  typedef CostVal REAL;
  DualDecompositionFlipperES(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecompositionFlipperES(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecompositionFlipperES();

  void setDDFlipperFunction(MRF::SmoothCostGeneralFn fnCost_dd_flipper);
  void setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper);
  void setNodeMappingStructure(std::vector<int> **sg_node_to_g_node);
  void setflippedinfo(std::vector<int> **flippedinfo);
  void set_edge_to_subgraphcount(myumap **edge_to_subgraphcount_p);


  // subnonsub breaks graph into one submodular graph and minimum number of nonsub subgraphs
  // for a planar grid graph, only one nonsubmod graph is sufficient
  // TREESADD starts with the tree splitting used by DD-T and then converts them to flipper graphs

  void setSplitType(DDFES::SplitType s_var);

protected:

  virtual void optimizeAlg(int nIterations);

  std::vector<int> **subgraph_node_to_graph_node_ptop;
  SmoothCostGeneralFn m_ddflippersmoothFn;
  SmoothCostGeneralFn m_flippersmoothFn;
  std::vector<int> **flippedinfo_ptop;
  myumap **edge_to_subgraphcount_ptop;

private:

  DDFES::SplitType m_splittype;

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

  void createSubgraphsTreeAddES(
      const std::vector<std::vector<std::vector<int> > > &node_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<int> > *subgraphs_edge_list,
      std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes);

  void addEdgesToFlipper(
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &edge_type,
      std::vector<int> *mark_flip_nodes,
      std::vector<int> *subgraph_edge_list);


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

  void updateCountEdgeSharing(
      const std::vector<std::vector<int> > &subgraphs_edge_list,
      std::vector<int> *edges_count_shared_subgraphs);

  void updateMapCountEdgeSharing(
      const std::vector<int> &edges_count_shared_subgraphs,
      const std::vector<std::vector<int> > &edge_list,
      myumap *edge_to_subgraphcount);

  void copy_marked_flipped_nodes(
      const std::vector<int> &mark_flip_nodes,
      std::vector<std::vector<int> > *mark_flip_nodes_sparse);


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



#endif /* DUALDECOMPOSITIONFLIPPERES_H_ */
