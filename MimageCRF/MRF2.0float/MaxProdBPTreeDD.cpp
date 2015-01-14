/*
 * MaxProdBPTreeDD.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: bhole
 */

#include "MaxProdBPTreeDD.h"


MaxProdBPTreeDD::MaxProdBPTreeDD(int nPixels, int nLabels,EnergyFunction *eng):MaxProdBPTree(nPixels, nLabels, eng)
{
  share_edges = 0;
  tree_node_to_graph_node = 0;
}

void MaxProdBPTreeDD::set_tree_node_to_graph_node(std::vector<int> *tngn)
{
  tree_node_to_graph_node = tngn;
}

void MaxProdBPTreeDD::set_edge_to_treecount(myumap *ettc)
{
  edge_to_treecount = ettc;
}

MRF::EnergyVal MaxProdBPTreeDD::call_smoothFn(int pix_i, int pix_j, int label_i, int label_j)
{
  // apply pixel transformation here.
  // convert pix1 and pix2 of the tree to pixels of the graph
  // say pixg_i and pixg_j
  int pixg_i = (*tree_node_to_graph_node)[pix_i];
  int pixg_j = (*tree_node_to_graph_node)[pix_j];

  MRF::SmoothCostGeneralFn function_name = MaxProdBPTree::getSmoothnessFn();

  if (share_edges == 0)
  {
    return function_name(pixg_i, pixg_j, label_i, label_j);//, graph_smoothFn, tree_node_to_graph_node);
  } else if (share_edges == 1) {
    int pix1, pix2; // pix1 smaller than pix2

    if (pixg_i <= pixg_j)
    {
      pix1 = pixg_i;
      pix2 = pixg_j;
    } else {
      pix2 = pixg_i;
      pix1 = pixg_j;
    }

    char buff[30];
    sprintf(buff, "%012d%012d", pix1, pix2);
    std::string str(buff);

    myumap::const_iterator got = (*edge_to_treecount).find (buff);

    int tree_count = 0;
    if ( got == (*edge_to_treecount).end() )
    {
      printf("\n Real problem here : pix1 %d pix2 %d  buff %s\n", pix1, pix2, buff);
      exit(1);
    }
    else
      tree_count = got->second;

    return double(function_name(pixg_i, pixg_j, label_i, label_j))/double(tree_count);

  }
}
