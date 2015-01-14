/*
 * DualDecompositionTree.h
 *
 *  Created on: Dec 28, 2011
 *      Author: bhole
 */

#ifndef DUALDECOMPOSITIONTREE_H_
#define DUALDECOMPOSITIONTREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
//#include "mrf.h"
#include "DualDecomposition.h"
//#include "LinkedBlockList.h"
//#include "regions-new.h"

#define FloatType float
#define FLOATTYPE float

//class OneNodeCluster;

class DualDecompositionTree : public DualDecomposition {
public:
  typedef CostVal REAL;
  DualDecompositionTree(int width, int height, int nLabels, EnergyFunction *eng);
  DualDecompositionTree(int nPixels, int nLabels,EnergyFunction *eng);
  ~DualDecompositionTree();

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


};



#endif /* DUALDECOMPOSITIONTREE_H_ */
