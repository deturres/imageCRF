/*
 * DualDecomposition.h
 *
 *  Created on: Oct 3, 2011
 *      Author: bhole
 */

#ifndef DUALDECOMPOSITION_H_
#define DUALDECOMPOSITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "MRFWrapper.h"


class DualDecomposition : public MRFWrapper {
public:
  typedef CostVal REAL;
  DualDecomposition(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecomposition(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecomposition();


  double getToleranceParam() { return primaldualtolerance; }
  void setToleranceParam(double tol) { primaldualtolerance = tol; }

protected:

  int getEdgeForNodes(const int &start_node,  const int &neighbour,
                      const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list);

  void extractDFSTreeMaxLimitNonRepeatableEdges(
      const int &max_edge_limit,
      const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *global_used_edges,
      std::vector<std::vector<int> > *trees_edge_list);

  int areThereUnusedEdges(const std::vector<int> &global_used_edges);

  void extractKruskalMaxLimitRepeatableEdges(
      const int &max_edge_limit,
      const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *global_used_edges,
      std::vector<std::vector<int> > *trees_edge_list);

  int attemptAddEdge(const int curr_edge,
      const std::vector<std::vector<int> > &graph_edge_list,
      std::vector<int> *nodes_in_what_tree,
      int *kruskal_forest_size,
      std::vector<std::vector <int> > *kruskal_forest,
      int *num_edges);



private:

  double primaldualtolerance;

};



#endif /* DUALDECOMPOSITION_H_ */
