/*
 * DualDecompositionTreeES.h
 *
 *  Created on: May 27, 2012
 *      Author: bhole
 */

// edge sharing enabled.

#ifndef DUALDECOMPOSITIONTREEES_H_
#define DUALDECOMPOSITIONTREEES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include "common_typedefs.h"
#include "DualDecomposition.h"

#define FloatType float
#define FLOATTYPE float


class DualDecompositionTreeES : public DualDecomposition {
public:
  typedef CostVal REAL;
  DualDecompositionTreeES(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecompositionTreeES(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecompositionTreeES();

protected:

  virtual void optimizeAlg(int nIterations);

private:

  void extractPrimTreeEdges(const int &prefer_not_used_edges,
                            const std::vector<std::vector<std::vector<int> > > &node_edge_list,
                            const std::vector<std::vector<int> > &edge_list,
                            std::vector<int> *global_used_edges,
                            std::vector<std::vector<int> > *trees_edge_list);

  void computeDataCostForTrees(
      const std::vector<std::vector<int> > &trees_edge_list,
      const std::vector<std::vector<int> > &edge_list,
      std::vector<std::vector<std::vector<int> > > &node_in_trees,
      std::vector<MRF::CostVal*> *D_a,
      std::vector<std::vector<int> > *tree_node_to_graph_node,
      std::vector<DataCost*> *dcost_vec);

  void updateCountEdgeSharing(
      const std::vector<std::vector<int> > &trees_edge_list,
      std::vector<int> *edges_count_shared_trees);

  void updateMapCountEdgeSharing(
      const std::vector<int> &edges_count_shared_trees,
      const std::vector<std::vector<int> > &edge_list,
      myumap *edge_to_treecount);

};



#endif /* DUALDECOMPOSITIONTREEES_H_ */
