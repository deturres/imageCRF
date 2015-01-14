/*
 * MaxProdBPTreeDD.h
 *
 *  Created on: Oct 15, 2011
 *      Author: bhole
 */

#ifndef MAXPRODBPTREEDD_H_
#define MAXPRODBPTREEDD_H_

#include "MaxProdBPTree.h"
#include "common_typedefs.h"

class MaxProdBPTreeDD : public MaxProdBPTree{

public:
  MaxProdBPTreeDD(int nPixels, int nLabels,EnergyFunction *eng);
  ~MaxProdBPTreeDD() {};
  virtual MRF::EnergyVal call_smoothFn(int pix_i, int pix_j, int label_i, int label_j);
  void set_tree_node_to_graph_node(std::vector<int> *tngn);
  void set_edge_to_treecount(myumap *ettc);
  void set_shared_edges(int value)
  { share_edges = value; }

private:
  myumap *edge_to_treecount;
  std::vector<int> *tree_node_to_graph_node;
  int share_edges;

};


#endif /* MAXPRODBPTREEDD_H_ */
