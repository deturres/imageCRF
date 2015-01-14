/*
 * DualDecompositionSubGraph.h
 *
 *  Created on: Dec 28, 2011
 *      Author: bhole
 */

#ifndef DUALDECOMPOSITIONSUBGRAPH_H_
#define DUALDECOMPOSITIONSUBGRAPH_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "DualDecomposition.h"


#define FloatType float
#define FLOATTYPE float

//class OneNodeCluster;

class DualDecompositionSubGraph : public DualDecomposition {
public:
  typedef CostVal REAL;
  DualDecompositionSubGraph(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecompositionSubGraph(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecompositionSubGraph();

  void setSubmodularFunction(MRF::SmoothCostGeneralFn subfnCost);
  void setNonSubmodularFunction(MRF::SmoothCostGeneralFn nonsubfnCost);
  void setNodeMappingStructure(std::vector<int> **sg_node_to_g_node);
  void setflippedinfo(std::vector<int> **flippedinfo);

protected:

  virtual void optimizeAlg(int nIterations);
  SmoothCostGeneralFn m_ssmoothFn; // submodular
  SmoothCostGeneralFn m_nssmoothFn; //non-submodular
  std::vector<int> **subgraph_node_to_graph_node_ptop;
  std::vector<int> **flippedinfo_ptop;

private:

  int checkSubmodularityEdge(const int &node1, const int &node2);

  void extractSubGraphEdges(const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *global_used_edges,
      std::vector<std::vector<int> > *subgraphs_edge_list,
      std::vector<int> *subgraph_type);

  void markFlippedNodes(const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &subgraph_type,
      std::vector<std::vector<std::vector<int> > >*node_in_subgraphs);

  void computeDataCostForSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      const std::vector<int> &subgraph_type,
      std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
      std::vector<MRF::CostVal*> *D_a,
      std::vector<std::vector<int> > *subgraph_node_to_graph_node,
      std::vector<DataCost*> *dcost_vec,
      std::vector<std::vector<int> > *flipstructure);

  void getSubGraphNodeNo(const std::vector<std::vector<int> > &node_in_subgraphs_g,
      const int &subgraph_no, int *subgraph_node1);

  void getWhatSubGraphsNodeIsUsedIn(
      const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
      std::vector<int> *num_nodes_in_subgraphs);

  void dumpSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
      const std::vector<std::vector<int> > &edge_list);


};


#endif /* DUALDECOMPOSITIONSUBGRAPH_H_ */
