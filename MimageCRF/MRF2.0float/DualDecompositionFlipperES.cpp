/*
 * DualDecompositionFlipperES.cpp
 *
 *  Created on: May 29, 2012
 *      Author: bhole
 */

// edge sharing enabled.


#include <time.h>
#include <limits.h>
#include <math.h>
#include <queue>

#include "DualDecompositionFlipperES.h"
#include "FlipperForDD.h"



DualDecompositionFlipperES::DualDecompositionFlipperES(int width, int height, int nLabels,EnergyFunction *eng):DualDecomposition(width,height,nLabels,eng)
{
}
DualDecompositionFlipperES::DualDecompositionFlipperES(int nPixels, int nLabels,EnergyFunction *eng):DualDecomposition(nPixels,nLabels,eng)
{
}

DualDecompositionFlipperES::~DualDecompositionFlipperES()
{
}

void DualDecompositionFlipperES::setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper)
{
  m_flippersmoothFn = fnCost_flipper;
}

void DualDecompositionFlipperES::setDDFlipperFunction(MRF::SmoothCostGeneralFn fnCost_dd_flipper)
{
  m_ddflippersmoothFn = fnCost_dd_flipper;
}


void DualDecompositionFlipperES::setflippedinfo(std::vector<int> **flippedinfo_f)
{
  flippedinfo_ptop = flippedinfo_f;
}


void DualDecompositionFlipperES::setNodeMappingStructure(
    std::vector<int> **sg_node_to_g_node)
{
  subgraph_node_to_graph_node_ptop = sg_node_to_g_node;
}


void DualDecompositionFlipperES::set_edge_to_subgraphcount(
    myumap **edge_to_subgraphcount_p)
{
  edge_to_subgraphcount_ptop = edge_to_subgraphcount_p;
}



void DualDecompositionFlipperES::getSubGraphNodeNo(const std::vector<std::vector<int> > &node_in_subgraphs_g,
                   const int &subgraph_no, int *subgraph_node1)
{
  *subgraph_node1 = -1;
  for (unsigned int i = 0; i < node_in_subgraphs_g.size(); i++)
  {
    if (subgraph_no == node_in_subgraphs_g[i][0])
    {
      *subgraph_node1 = node_in_subgraphs_g[i][1];
    }
  }
  assert(*subgraph_node1 != -1);
}



void DualDecompositionFlipperES::setSplitType(DDFES::SplitType s_var)
{
  m_splittype = s_var;
}




void DualDecompositionFlipperES::createSubgraphs(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector<int> > *subgraphs_edge_list,
    std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes)
{
  if (m_splittype == DDFES::SUBNONSUB)
  {

  } else if (m_splittype == DDFES::TREESADD)
  {
    createSubgraphsTreeAdd(node_edge_list, edge_list, subgraphs_edge_list, subgraphs_flipped_nodes);

  } else if (m_splittype == DDFES::TREESADD_EDGESHARE)
  {
    createSubgraphsTreeAddES(node_edge_list, edge_list, subgraphs_edge_list, subgraphs_flipped_nodes);
  }
}



// assume for the spanning tree ie first subgraph, no need for checking all nodes in an external
// loop that is used for forest of graphs
void DualDecompositionFlipperES::markFlippedNodesFlipper(
    const std::vector<int> &subgraph_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &edge_type,
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    std::vector<int> *mark_flip_nodes)
{
  if (subgraph_edge_list.size() <= 0)
  {
      printf(" Subgraph empty, something fishy, lobstery, prawny!! ");
      exit(1);
  }

  if (subgraph_edge_list.size() > 0) // spanning tree has been copied to subgraph
  {

    // so that we can have a O(1) call to find out what edges exist in current flipper graph
    std::vector<int> current_edge_list(edge_list.size(), 0);

    for (int i=0; i<subgraph_edge_list.size(); i++)
    {
      int e = subgraph_edge_list[i];
      current_edge_list[e] = 1;
    }

    std::queue<int> flippedlist;
    std::vector<int> checked(m_nPixels, 0);

    int firstedge = subgraph_edge_list[0];
    int firsttype = edge_type[firstedge];
    int node1 = edge_list[firstedge][0];
    int node2 = edge_list[firstedge][1];
    if (firsttype == 0) // submod
    {
      (*mark_flip_nodes)[node1] = 1; // not flipped
      (*mark_flip_nodes)[node2] = 1; // not flipped
    } else if (firsttype == 1) // nonsubmod
    {
      (*mark_flip_nodes)[node1] = 1; // not flipped
      (*mark_flip_nodes)[node2] = 2; // flipped
    }
    flippedlist.push(node1);
    flippedlist.push(node2);

    checked[node1] = 1;
    checked[node2] = 1;

    // printf(" Initial nodes : %d, %d , first edge : %d \n", node1, node2, firstedge);

    while (!flippedlist.empty())
    {
      // nodep stands for parent node
      int nodep = flippedlist.front();
      flippedlist.pop();

      for (int nn=0; nn < node_edge_list[nodep].size(); nn++)
      { // go through all neighbors
        int neighbor = node_edge_list[nodep][nn][0];
        int edge = node_edge_list[nodep][nn][1];

        if (current_edge_list[edge] == 1) // this exists in the current tree/graph
        {
          if (checked[neighbor] == 0) // not visited yet, so can add this edge
          {
            int ftype = edge_type[edge];
            if (ftype == 0) // submod - both nodes are either flipped or not flipped
            {
              if ((*mark_flip_nodes)[nodep] == 1)
                (*mark_flip_nodes)[neighbor] = 1; // not flipped

              if ((*mark_flip_nodes)[nodep] == 2)
                (*mark_flip_nodes)[neighbor] = 2; // flipped

            } else if (ftype == 1) // nonsubmod - one node is flipped, other not flipped
            {
              if ((*mark_flip_nodes)[nodep] == 1)
                (*mark_flip_nodes)[neighbor] = 2; // flipped

              if ((*mark_flip_nodes)[nodep] == 2)
                (*mark_flip_nodes)[neighbor] = 1; // not flipped
            }

            flippedlist.push(neighbor);
            checked[neighbor] = 1;

            // printf("Parent node : %d => Neighbor :%d, Edge : %d \n", nodep, neighbor, edge);

          }
        }
      }
    }
  }
}




void DualDecompositionFlipperES::moveTreeToSubgraph(
    std::vector<int> *tree_edge_list,
    std::vector<std::vector<int> > *subgraphs_edge_list)
{

  std::vector<int> subgraph_edge_list;

  for (int e=0; e <tree_edge_list->size(); e++)
  {
    subgraph_edge_list.push_back((*tree_edge_list)[e]);
  }

  subgraphs_edge_list->push_back(subgraph_edge_list);

  // erase tree edges
  tree_edge_list->erase(tree_edge_list->begin(), tree_edge_list->end());

}


void DualDecompositionFlipperES::traverseGraph(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<int> &subgraphs_edge_list0,
    const std::vector<std::vector<int> > &edge_list)
{

  if (subgraphs_edge_list0.size() <= 0 )
  {
    printf(" Empty graph \n");
    return;
  }

  std::vector<int> current_edge_list(edge_list.size(), 0);

  for (int i=0; i<subgraphs_edge_list0.size(); i++)
  {
    int e = subgraphs_edge_list0[i];
    current_edge_list[e] = 1;
  }

  std::queue<int> traverselist;
  std::vector<int> checked(m_nPixels, 0);
  std::vector<int> shown(edge_list.size(), 0); // edges that are already displayed

  int firstedge = subgraphs_edge_list0[0];
  int node1 = edge_list[firstedge][0];
  int node2 = edge_list[firstedge][1];

  checked[node1] = 1;
  checked[node2] = 1;

  traverselist.push(node1);
  traverselist.push(node2);

  // printf(" Initial nodes : %d, %d , first edge : %d \n", node1, node2, firstedge);
  printf("Parent node : %d => Neighbor :%d, Edge : %d \n", node1, node2, firstedge);

  shown[firstedge] = 1;

  int alldone = 0;

  while (alldone == 0)
  {
    while (!traverselist.empty())
    {
      // nodep stands for parent node
      int nodep = traverselist.front();
      traverselist.pop();

      for (int nn=0; nn < node_edge_list[nodep].size(); nn++)
      { // go through all neighbors
        int neighbor = node_edge_list[nodep][nn][0];
        int edge = node_edge_list[nodep][nn][1];

        if (current_edge_list[edge] == 1) // this exists in the current tree/graph
        {
          if (checked[neighbor] == 0) // not visited yet, so can add node
          {
            traverselist.push(neighbor);
            checked[neighbor] = 1;

            printf("Parent node : %d => Neighbor :%d, Edge : %d \n", nodep, neighbor, edge);
            shown[edge] = 1;
          } else if (checked[neighbor] == 1 && shown[edge] == 0) // visited, only show edge
          {
            printf("Parent node : %d => Neighbor :%d, Edge : %d \n", nodep, neighbor, edge);
            shown[edge] = 1;
          }
        }

      }
    }

    alldone = 1;
    for (int k=0; k<current_edge_list.size(); k++)
    {
      if (current_edge_list[k] == 1 && shown[k]==0)
      {
        int node1 = edge_list[k][0];
        int node2 = edge_list[k][1];

        checked[node1] = 1;
        checked[node2] = 1;

        traverselist.push(node1);
        traverselist.push(node2);

        // printf(" Initial nodes : %d, %d , first edge : %d \n", node1, node2, k);
        printf("Parent node : %d => Neighbor :%d, Edge : %d \n", node1, node2, k);

        shown[k] = 1;
        alldone = 0;
        break;
      }
    }
  }
}


void DualDecompositionFlipperES::dumpSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
               const std::vector<std::vector<int> > &edge_list)
{
  for (unsigned int i=0; i < subgraphs_edge_list.size(); i++)
  {
    printf("\n SubGraph : %d SubGraph Edge size : %d\n", i, subgraphs_edge_list[i].size());
    for (unsigned int j = 0; j < subgraphs_edge_list[i].size(); j++)
    {
      int edge = subgraphs_edge_list[i][j];
      if (edge == -1)
      {
        printf(" Error in edge ");
        exit(1);
      }
      printf("Edge: %d Node: %d Node: %d\n", edge, edge_list[edge][0], edge_list[edge][1]);
    }
  }
}




void DualDecompositionFlipperES::displayNodeEdgeType(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &edge_type)
{

  for (int i=0; i < node_edge_list.size(); i++)
  {
    for (int j=0; j<node_edge_list[i].size(); j++)
    {
      printf(" Node : %d Neighbor : %d, Edge : %d Type : %d\n", i, node_edge_list[i][j][0],
          node_edge_list[i][j][1] , edge_type[node_edge_list[i][j][1]] );
    }
  }


}

// for now we try to insert only in subgraphs0
void DualDecompositionFlipperES::moveTreeEdgesToFlipper(
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &edge_type,
    std::vector<int> *mark_flip_nodes,
    std::vector<int> *subgraph_edge_list,
    std::vector<int> *tree_edge_list)
{

  std::vector<int> to_delete_indices;

  for (int i=0; i<tree_edge_list->size(); i++)
  {
    int edge = (*tree_edge_list)[i];
    // attempt to insert edge into flipper graph
    int etype = edge_type[edge];

    int node1 = edge_list[edge][0];
    int node2 = edge_list[edge][1];

    // printf(" Tree: edge %d edge type %d node1 %d node2 %d \n", edge, etype, node1, node2);
    // printf(" mark_node1 %d mark_node2 %d\n", (*mark_flip_nodes)[node1], (*mark_flip_nodes)[node2]);

    if ((*mark_flip_nodes)[node1] == 0 && (*mark_flip_nodes)[node2] == 0)
    {
      // both nodes not in flipper graph, add both of them and edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);

      if (etype == 0) // submod
      {
        (*mark_flip_nodes)[node1] = 1;
        (*mark_flip_nodes)[node2] = 1;
      } else if (etype == 1) // nonsubmod
      {
        (*mark_flip_nodes)[node1] = 1;
        (*mark_flip_nodes)[node2] = 2;
      }

    } else if ((*mark_flip_nodes)[node1] == 0 && (*mark_flip_nodes)[node2] == 1) // node2 is not flipped
    {
      // add node1 and edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
      if (etype == 0) // submod
      {
        (*mark_flip_nodes)[node1] = 1;
      } else if (etype == 1) // nonsubmod
      {
        (*mark_flip_nodes)[node1] = 2;
      }

    } else if ((*mark_flip_nodes)[node1] == 0 && (*mark_flip_nodes)[node2] == 2) // node2 is flipped
    {
      // add node1 and edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
      if (etype == 0) // submod
      {
        (*mark_flip_nodes)[node1] = 2;
      } else if (etype == 1) // nonsubmod
      {
        (*mark_flip_nodes)[node1] = 1;
      }
    } else if ((*mark_flip_nodes)[node2] == 0 && (*mark_flip_nodes)[node1] == 1)
    {
      // add node2 and edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
      if (etype == 0) // submod
      {
        (*mark_flip_nodes)[node2] = 1;
      } else if (etype == 1) // nonsubmod
      {
        (*mark_flip_nodes)[node2] = 2;
      }
    } else if ((*mark_flip_nodes)[node2] == 0 && (*mark_flip_nodes)[node1] == 2)
    {
      // add node2 and edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
      if (etype == 0) // submod
      {
        (*mark_flip_nodes)[node2] = 2;
      } else if (etype == 1) // nonsubmod
      {
        (*mark_flip_nodes)[node2] = 1;
      }
    } else if ((*mark_flip_nodes)[node1] == 1 && (*mark_flip_nodes)[node2] == 1 && etype == 0)
    {
      // add edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
    } else if ((*mark_flip_nodes)[node1] == 2 && (*mark_flip_nodes)[node2] == 2 && etype == 0)
    {
      // add edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
    } else if ((*mark_flip_nodes)[node1] == 1 && (*mark_flip_nodes)[node2] == 2 && etype == 1)
    {
      // add edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
    } else if ((*mark_flip_nodes)[node1] == 2 && (*mark_flip_nodes)[node2] == 1 && etype == 1)
    {
      // add edge
      subgraph_edge_list->push_back(edge);
      to_delete_indices.push_back(i);
    }

  }

  // erase from largest index to maintain correctness of erasing vector
  for (int i=to_delete_indices.size()-1; i>=0; i--)
  {
    tree_edge_list->erase(tree_edge_list->begin()+to_delete_indices[i]);
  }

}




// for now we try to insert only in subgraphs0
void DualDecompositionFlipperES::addEdgesToFlipper(
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &edge_type,
    std::vector<int> *mark_flip_nodes,
    std::vector<int> *subgraph_edge_list)
{
  // this stores whether edge is contained in subgraph or not, it is a dense substitute for the
  // sparse subgraph_edge_list
  std::vector<int> subgraph_edge_list_array(edge_list.size(), 0);

  for (int i=0; i<(*subgraph_edge_list).size(); i++)
  {
    subgraph_edge_list_array[(*subgraph_edge_list)[i]] = 1;
  }


  for (int e=0; e<edge_list.size(); e++)
  {
    // make sure that edge e is already not in the subgraph
    if (subgraph_edge_list_array[e] == 1) // already exists
    {}
    else
    {
      int etype = edge_type[e];

      int node1 = edge_list[e][0];
      int node2 = edge_list[e][1];

      if ((*mark_flip_nodes)[node1] == 0 || (*mark_flip_nodes)[node2] == 0)
      {
        // since we start with spanning tree, all nodes must exist in the subgraph
        printf(" Nodes have to exist in the subgraph from before \n");
        exit(1);
      } else {

        if ((*mark_flip_nodes)[node1] == 1 && (*mark_flip_nodes)[node2] == 1 && etype == 0)
        {
          // add edge
          subgraph_edge_list->push_back(e);
        } else if ((*mark_flip_nodes)[node1] == 2 && (*mark_flip_nodes)[node2] == 2 && etype == 0)
        {
          // add edge
          subgraph_edge_list->push_back(e);
        } else if ((*mark_flip_nodes)[node1] == 1 && (*mark_flip_nodes)[node2] == 2 && etype == 1)
        {
          // add edge
          subgraph_edge_list->push_back(e);
        } else if ((*mark_flip_nodes)[node1] == 2 && (*mark_flip_nodes)[node2] == 1 && etype == 1)
        {
          // add edge
          subgraph_edge_list->push_back(e);
        }

      }
    }
  }

}











void DualDecompositionFlipperES::copy_marked_flipped_nodes(
    const std::vector<int> &mark_flip_nodes,
    std::vector<std::vector<int> > *mark_flip_nodes_sparse)
{
  for (int i=0; i < mark_flip_nodes.size(); i++)
  {
    if (mark_flip_nodes[i]!=0) // ie so it is used in flipper graph
    {
      std::vector<int> temp_struct(2, 0);
      temp_struct[0]=i; // storing node number
      if (mark_flip_nodes[i]==1) // not flipped
      {
        temp_struct[1]=1;
      }
      else if (mark_flip_nodes[i]==2) // flipped
      {
        temp_struct[1]=2;
      }
      mark_flip_nodes_sparse->push_back(temp_struct);
    }
  }
}

void DualDecompositionFlipperES::createSubgraphsTreeAdd(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector<int> > *subgraphs_edge_list,
    std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes
)
{

  // 0 indicates not used, 1 indicates edge is used by at least one tree.
  std::vector<int> global_used_edges(edge_list.size(), 0);

  // [number of trees][edges numbers]
  // ie [tree_no [edge_no]*]*
  std::vector<std::vector<int> > trees_edge_list;

  printf(" Creating Intial Basic Trees ... \n");
  while (areThereUnusedEdges(global_used_edges))
  {
    extractDFSTreeMaxLimitNonRepeatableEdges(edge_list.size()+1, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);
  }

  std::vector<int> edge_type(edge_list.size(), 0); // 0 stands for submod, 1 for nonsubmod
  extractEdgeType(edge_list, &edge_type);

  //displayNodeEdgeType(node_edge_list, edge_list, edge_type);

  int numtrees = trees_edge_list.size();

  int j = 0;
  while (j < numtrees)
  {

    // move tree edges to subgraph
    moveTreeToSubgraph(&trees_edge_list[j], subgraphs_edge_list);

    // 0 in this array indicates node not used in current flipper graph,
    // 1 indicates node is used in current flipper graph but is not flipped
    // 2 indicates node is used in current flipper graph and is flipped
    std::vector<int> mark_flip_nodes(m_nPixels, 0);

    // run bfs on all trees in current flipper graph (which could be a forest)
    // to mark used , unflipped and flipped nodes.
    int curr_flipper_graph = subgraphs_edge_list->size() - 1;
    // printf("current : %d\n", curr_flipper_graph);
    markFlippedNodesFlipper((*subgraphs_edge_list)[curr_flipper_graph], edge_list, edge_type, node_edge_list, &mark_flip_nodes);


    for (int k=j+1; k<numtrees; k++)
    {
      moveTreeEdgesToFlipper(edge_list, edge_type, &mark_flip_nodes, &(*subgraphs_edge_list)[curr_flipper_graph], &trees_edge_list[k]);
    }

    // copy marked flipped nodes so we don't need to do bfs again.
    std::vector<std::vector<int> > mark_flip_nodes_sparse; // [node_no flip]
    copy_marked_flipped_nodes(mark_flip_nodes, &mark_flip_nodes_sparse);
    subgraphs_flipped_nodes->push_back(mark_flip_nodes_sparse);

    j++;

    while (j < numtrees && trees_edge_list[j].size() == 0)
      j++;
  }

//  for (int k=0; k < subgraphs_edge_list->size(); k++)
//  {
//    printf(" Subgraph no. %d \n", k);
//    traverseGraph(node_edge_list, (*subgraphs_edge_list)[k], edge_list);
//  }


}





void DualDecompositionFlipperES::createSubgraphsTreeAddES(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector<int> > *subgraphs_edge_list,
    std::vector<std::vector<std::vector<int> > > *subgraphs_flipped_nodes
)
{

  // 0 indicates not used, 1 indicates edge is used by at least one tree.
  std::vector<int> global_used_edges(edge_list.size(), 0);

  // [number of trees][edges numbers]
  // ie [tree_no [edge_no]*]*
  std::vector<std::vector<int> > trees_edge_list;

  printf(" Creating Intial Basic Trees ... \n");

  extractDFSTreeMaxLimitNonRepeatableEdges(edge_list.size()+1, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);

  while (areThereUnusedEdges(global_used_edges))
  {
    extractKruskalMaxLimitRepeatableEdges(edge_list.size()+1, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);
  }
  printf(" Number of trees created %d \n", trees_edge_list.size());



  std::vector<int> edge_type(edge_list.size(), 0); // 0 stands for submod, 1 for nonsubmod
  extractEdgeType(edge_list, &edge_type);

  // displayNodeEdgeType(node_edge_list, edge_list, edge_type);

  int numtrees = trees_edge_list.size();

  int j = 0;
  while (j < numtrees)
  {
    // move spannning tree edges to subgraph
    moveTreeToSubgraph(&trees_edge_list[j], subgraphs_edge_list);
    j++;
  }


  // Attempt to add edges of graph in all the subgraphs so as to maintain flipper properties
  j = 0;
  while (j < numtrees)
  {

    // 0 in this array indicates node not used in current flipper graph => not possible since we start with spanning tree,
    // 1 indicates node is used in current flipper graph but is not flipped
    // 2 indicates node is used in current flipper graph and is flipped
    std::vector<int> mark_flip_nodes(m_nPixels, 0);

    // run bfs on all trees in current flipper graph (which could be a forest)
    // to mark used , unflipped and flipped nodes.
    int curr_flipper_graph = j;
    // printf("current : %d\n", curr_flipper_graph);
    markFlippedNodesFlipper((*subgraphs_edge_list)[curr_flipper_graph], edge_list, edge_type, node_edge_list, &mark_flip_nodes);

    addEdgesToFlipper(edge_list, edge_type, &mark_flip_nodes, &(*subgraphs_edge_list)[curr_flipper_graph]);

    // copy marked flipped nodes so we don't need to do bfs again.
    std::vector<std::vector<int> > mark_flip_nodes_sparse; // [node_no flip]
    copy_marked_flipped_nodes(mark_flip_nodes, &mark_flip_nodes_sparse);
    subgraphs_flipped_nodes->push_back(mark_flip_nodes_sparse);

    j++;
  }
/*
  for (int k=0; k < subgraphs_edge_list->size(); k++)
  {
    printf(" Subgraph no. %d \n", k);
    traverseGraph(node_edge_list, (*subgraphs_edge_list)[k], edge_list);
  }
*/

}







void DualDecompositionFlipperES::getWhatSubGraphsNodeIsUsedIn(
    const std::vector<std::vector<int> > &subgraphs_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
    std::vector<int> *num_nodes_in_subgraphs)
{
  const int num_subgraphs = subgraphs_edge_list.size();
  int edge_number;

  for (int t=0; t < num_subgraphs; t++)
  {
    for (unsigned int e=0; e < subgraphs_edge_list[t].size(); e++)
    {
      edge_number = subgraphs_edge_list[t][e];
      int node1 = edge_list[edge_number][0];
      int node2 = edge_list[edge_number][1];

      // for node1
      int subgraphs_size = (*node_in_subgraphs)[node1].size();
      int found = 0;
      for (int f=0; f<subgraphs_size; f++)
      {
        if (t == (*node_in_subgraphs)[node1][f][0])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
      {
        std::vector<int> inn(3);
        inn[0] = t;
        inn[1] = -1; // subgraph node number will be filled in later function
        inn[2] = 0; // default indicates not flipped
        (*node_in_subgraphs)[node1].push_back(inn);
        (*num_nodes_in_subgraphs)[t]++;
      }

      // for node2
      subgraphs_size = (*node_in_subgraphs)[node2].size();
      found = 0;
      for (int f=0; f<subgraphs_size; f++)
      {
        if (t == (*node_in_subgraphs)[node2][f][0])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
      {
        std::vector<int> inn(3);
        inn[0] = t;
        inn[1] = -1; // subgraph node number will be filled in later function
        inn[2] = 0; // default indicates not flipped
        (*node_in_subgraphs)[node2].push_back(inn);
        (*num_nodes_in_subgraphs)[t]++;
      }
    }
  }
}



/*

void updateNodeDataStructure(
    const std::vector<std::vector<std::vector<int> > > &subgraphs_flipped_nodes,
    std::vector<std::vector<std::vector<int> > > *node_in_subgraphs)
{
  // for each graph
  for (int g=0; g<subgraphs_flipped_nodes.size(); g++)
  {
    // for each node in the subgraph
    for (int i=0; i<subgraphs_flipped_nodes[g].size(); i++)
    {
      int nodem = subgraphs_flipped_nodes[g][i][0];
      int flipinfo = subgraphs_flipped_nodes[g][i][1];
      // flipinfo is 1 if not flipped and 2 if flipped

      // find subgraph in (*node_in_subgraphs)[nodem]
      for (int kk=0; kk<(*node_in_subgraphs)[nodem].size(); kk++)
      {
        if ((*node_in_subgraphs)[nodem][kk][0] == g)
        {
          if (flipinfo == 1)
          (*node_in_subgraphs)[nodem][kk][2] = 0; // not flipped
          else if (flipinfo == 2)
            (*node_in_subgraphs)[nodem][kk][2] = 1; // flipped
        }
      }
    }
  }
}

*/




// the node number order is maintained.
void DualDecompositionFlipperES::computeDataCostForSubGraphs(
    const std::vector<std::vector<int> > &subgraphs_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
    std::vector<MRF::CostVal*> *D_a,
    std::vector<std::vector<int> > *subgraph_node_to_graph_node,
    std::vector<DataCost*> *dcost_vec)
{

  assert(m_nLabels == 2);

  const int num_subgraphs = subgraphs_edge_list.size();

  std::vector<int> dsiIndex(num_subgraphs, 0);

  // node_in_subgraphs gives all the subgraphs that a node belongs to
  for (unsigned int i=0; i < node_in_subgraphs->size(); i++)
  {
    assert((*node_in_subgraphs)[i].size() == num_subgraphs);
    for (unsigned int t = 0; t < (*node_in_subgraphs)[i].size(); t++)
    {
      int subgraph_num = (*node_in_subgraphs)[i][t][0];
      (*subgraph_node_to_graph_node)[subgraph_num].push_back(i);
      int subgraph_node_no = (*subgraph_node_to_graph_node)[subgraph_num].size() - 1;
      (*node_in_subgraphs)[i][t][1] = subgraph_node_no;

      for (int d = 0; d < m_nLabels; d++)
      {
        MRF::CostVal* ptr = (*D_a)[subgraph_num] + dsiIndex[subgraph_num];
        *ptr = nodeArray[i].localEv[d] / (*node_in_subgraphs)[i].size();
        dsiIndex[subgraph_num]++;
      }
    }
  }

  for (unsigned int t = 0; t < subgraphs_edge_list.size(); t++)
  {
    (*dcost_vec)[t] = new DataCost ((*D_a)[t]);
  }
}


// not necessary to store which subgraphs but rather just the count
// because we will use that for fnCost/num_shared_trees
void DualDecompositionFlipperES::updateCountEdgeSharing(
    const std::vector<std::vector<int> > &subgraphs_edge_list,
    std::vector<int> *edges_count_shared_subgraphs)
{
  for (int t=0; t<subgraphs_edge_list.size(); t++)
  {
    for (int e=0; e<subgraphs_edge_list[t].size(); e++)
    {
      int edge = subgraphs_edge_list[t][e];
      (*edges_count_shared_subgraphs)[edge] = (*edges_count_shared_subgraphs)[edge] + 1;
    }
  }
}




void DualDecompositionFlipperES::updateMapCountEdgeSharing(
    const std::vector<int> &edges_count_shared_subgraphs,
    const std::vector<std::vector<int> > &edge_list,
    myumap *edge_to_subgraphcount)
{

  for (int e=0; e<edges_count_shared_subgraphs.size(); e++)
  {
    int node1 = edge_list[e][0];
    int node2 = edge_list[e][1];

    // pix1 will contain the smaller pixelno
    int pix1, pix2;

    if (node1 <= node2)
    {
      pix1 = node1;
      pix2 = node2;
    } else {
      pix1 = node2;
      pix2 = node1;
    }

    char buff[30];
    sprintf(buff, "%012d%012d", pix1, pix2);
    std::string str(buff);

    int num_shared_subgraphs = edges_count_shared_subgraphs[e];

    std::pair<std::string, int> pairsi(buff, num_shared_subgraphs);
    (*edge_to_subgraphcount).insert(pairsi);
  }

}





void DualDecompositionFlipperES::optimizeAlg(int nIterations)
{
  if (m_grid_graph)
  {
    printf(" Dual Decomposition Flipper ES does not handle grid graphs, use other interface \n");
    exit(1);
  } else {
    // for non-grid graphs

    // Until number of iterations or convergence.
    // 1. use current graph to create subgraphs
    // 2. Split datacost and smoothcost terms for subgraphs
    // 3. Set up flipper pointers
    // 4. use Flipper code for each subgraph (after flipping trick) to do inference
    // 5. update parameters using projected sub-gradients.

    // srand ( time(NULL) );
    //srand(0);

    // [node_no, [other_no, edge_no]*]*
    std::vector<std::vector<std::vector<int> > > node_edge_list(m_nPixels);

    // [edge_no, [first_node, second_node]]*
    std::vector<std::vector<int> > edge_list;

    extractEdgeList(&node_edge_list, &edge_list);

    // 0 indicates not used, 1 indicates edge is used by at least one tree.
    std::vector<int> global_used_edges(edge_list.size(), 0);

    // ie [subgraph_no [edge_no]*]*
    std::vector<std::vector<int> > subgraphs_edge_list;

    if (m_nLabels != 2)
    {
      printf(" Current algorithm does not work for other than binary case. \n");
      return;
    }

    printf(" Creating Subgraphs ... \n");
    // indexed by subgraph no, then each node, not flipped(value=1) or flipped (value=2)
    // [subgraph no [node_no flipp]*]*
    std::vector<std::vector<std::vector<int> > > subgraphs_flipped_nodes;
    createSubgraphs(node_edge_list, edge_list, &subgraphs_edge_list, &subgraphs_flipped_nodes);

    // subgraphs_flipped_nodes will only provide proof as what nodes are flipped and only for debugging
    // purposes, otherwise these values are not needed as Flipper() call will handle everything
    // like a black box

    // after calling all subgraphs dump edges for each subgraph for debugging purpose.
    // dumpSubGraphs(subgraphs_edge_list, edge_list);

    const int num_subgraphs = subgraphs_edge_list.size();


    // this data structure stores how many trees a particular edge is shared among
    std::vector<int> edges_count_shared_subgraphs(edge_list.size(), 0);
    updateCountEdgeSharing(subgraphs_edge_list, &edges_count_shared_subgraphs);

    // edge is stored as a string [pixel1pixel2] with pixel1 being smaller
    // this is so that modified version of fnCost can find it easy to find the treecount value
    myumap edge_to_subgraphcount;
    updateMapCountEdgeSharing(edges_count_shared_subgraphs, edge_list, &edge_to_subgraphcount);
    *edge_to_subgraphcount_ptop = &edge_to_subgraphcount;

    // It is important to note the node numbering used for subgraphs.
    // ie what node numbers in the subgraph correspond to what node numbers
    // of the full graph
    // same for the edges.
    // hence we have mapping functions for both these.
    // the node number mapping is got from computeDataCostForTree

    // flip variable is used for indicating what node in the subgraph is flipped
    // this data structure is [graph_node_no][subgraph][0:subgraph_no, 1:subgraph_node_no 2:flipped]
    std::vector<std::vector<std::vector<int> > > node_in_subgraphs(m_nPixels, (std::vector<std::vector<int> >)0);

    // this structure counts number of nodes in each subgraph
    std::vector<int> num_nodes_in_subgraphs(num_subgraphs, 0);

    getWhatSubGraphsNodeIsUsedIn(subgraphs_edge_list, edge_list, &node_in_subgraphs, &num_nodes_in_subgraphs);

    for (unsigned int i=0; i < node_in_subgraphs.size(); i++)
    {
      // make sure that each node is used in all subgraphs.
      assert(node_in_subgraphs[i].size() == subgraphs_edge_list.size());
    }


    std::vector<DataCost*> dcost_vec(num_subgraphs, (DataCost*)0);

    std::vector<EnergyFunction*> energy_vec(num_subgraphs, (EnergyFunction*)0);
    std::vector<MRF*> mrf_vec(num_subgraphs, (MRF*)0);

    std::vector<MRF::CostVal*> D_a(subgraphs_edge_list.size(), (MRF::CostVal*)0);
    for (unsigned int t=0; t<num_subgraphs; t++)
    {
      D_a[t] = new MRF::CostVal[(num_nodes_in_subgraphs[t])*m_nLabels];
    }

    std::vector<std::vector<int> > subgraph_node_to_graph_node(num_subgraphs, (std::vector<int>)0);

    // creating flipping datastructure because it will be easier to use for flippinginfo
    // note that the 2nd index is on the subgraph node number, so that vector size varies
    // NOT NEEDED
    // std::vector<std::vector<int> > flipstructure(num_subgraphs, (std::vector<int>)0);

    // node_in_subgraphs contains a field with flip information
    // which needs to be updated
    // NOT NEEDED - so now node_in_subgraphs with flipper will contain invalid values
    // updateNodeDataStructure(subgraphs_flipped_nodes, &node_in_subgraphs);


    printf(" Computing DataCost for subgraphs ... \n");
    computeDataCostForSubGraphs(subgraphs_edge_list, edge_list, &node_in_subgraphs, &D_a, &subgraph_node_to_graph_node, &dcost_vec);
    // datacost with flipped nodes have datacost terms flipped.
    // so we need to remember that flipping while changing the parameters

    //this needs to call ddflipper because we need to do mapping
    SmoothnessCost* scost_common = new SmoothnessCost((MRF::SmoothCostGeneralFn) m_ddflippersmoothFn); // main function

    // create mrfs of Expansion for each subgraph we have.
    for (int t=0; t < num_subgraphs; t++)
    {
      energy_vec[t] = new EnergyFunction(dcost_vec[t], scost_common);

      // printf(" Size %d \n", subgraph_node_to_graph_node[t].size());
      mrf_vec[t] = new FlipperForDD(subgraph_node_to_graph_node[t].size(), m_nLabels, energy_vec[t]);

      // you need to set submodular/nonsubmodular function
      ((FlipperForDD*)mrf_vec[t])->setFlipperFunction((MRF::SmoothCostGeneralFn) m_flippersmoothFn);

      // set the edges for mrfs now
      for (int e=0; e < subgraphs_edge_list[t].size(); e++)
      {
        int graph_edge_no = subgraphs_edge_list[t][e];
        int graph_node1 = edge_list[graph_edge_no][0];
        int graph_node2 = edge_list[graph_edge_no][1];
        int subgraph_node1, subgraph_node2;
        getSubGraphNodeNo(node_in_subgraphs[graph_node1], t, &subgraph_node1);
        getSubGraphNodeNo(node_in_subgraphs[graph_node2], t, &subgraph_node2);
        mrf_vec[t]->setNeighbors(subgraph_node1, subgraph_node2, 1);
      }

      // Check validity of nodes connected in mrf.
      // will need to work on this.

      mrf_vec[t]->initialize();
      mrf_vec[t]->clearAnswer();

      ((FlipperForDD*)mrf_vec[t])->setflippedinfo(flippedinfo_ptop);
      *subgraph_node_to_graph_node_ptop = &subgraph_node_to_graph_node[t];

      ((FlipperForDD*)mrf_vec[t])->setUpFlipperForDD();

    }

    // stop timer
    clock_t finishss = clock();
    float timess = (float) (((double)(finishss - start)) / CLOCKS_PER_SEC);
    printf(" Time for creating and initializing subgraphs = %f secs\n", timess);

    // VERY IMPORTANT
    // Every time you call mrf_vec[t]->optimize(1, time_t); or mrf_vec[t]->totalEnergy();
    // you need to make sure that
    // the correct flipstructure and subgraph_node_to_graph_node structure are activated
    // for the external function


    int changed = 1;
    int iter_no = 1;
    int max_iter = 400; // refers to infer from subgraph not gradient descent times
    // ie gradient descent times is 1 less than max_iter

    float gamma = 0.1;

    int gmismatchno = m_nPixels * num_subgraphs;
    std::vector<int> glabels(m_nPixels);


    MRF::EnergyVal e_subgraphs_total_old = 0.0;
    MRF::EnergyVal e_subgraphs_total = 0.0;
    MRF::EnergyVal best_primal_energy = 0.0;

    printf("\n*******  Iterating ****\n");



    while (changed == 1)
    {
      changed = 0;
      // when any label in any tree changes, set changed to 1.

      e_subgraphs_total = 0.0;



      for (int t=0; t < num_subgraphs; t++)
      {
        MRF::EnergyVal E;
        float time_t, tot_time_t;
        // *flippedinfo_ptop = &flipstructure[t];
        // call before calling calls to Energy or optimize for each graph
        // you need to set *flippedinfo before calling that function
        ((FlipperForDD*)mrf_vec[t])->setflippedinfo(flippedinfo_ptop); //&flippedinfo_f);
        *subgraph_node_to_graph_node_ptop = &subgraph_node_to_graph_node[t];
/*
        for (int kk=0; kk<(**subgraph_node_to_graph_node_ptop).size(); kk++)
        {
          printf("subgraph node no %d graph node no %d\n", kk, (**subgraph_node_to_graph_node_ptop)[kk]);
        }

        for (int kk=0; kk<subgraph_node_to_graph_node[t].size(); kk++)
        {
          printf("subgraph node no %d graph node no %d\n", kk, subgraph_node_to_graph_node[t][kk]);
        }
*/

        E = mrf_vec[t]->totalEnergy();
        tot_time_t= 0;
        mrf_vec[t]->optimize(1, time_t);
        E = mrf_vec[t]->totalEnergy();
        tot_time_t = tot_time_t + time_t ;

//        for (int ix=0; ix<num_nodes_in_subgraphs[t]; ix++)
//        {
//          unsigned char rowvalue = mrf_vec[t]->getLabel(ix);
//          printf("Tree : %d  Node : %d  state : %d \n", t, ix, (int)rowvalue);
//        }


        e_subgraphs_total += E;
      }

      // printf("**** Computing norm results from the trees ****\n");
      // I call the sum of result labels divide by number of trees that contain the node
      // as norm_results.
      std::vector<std::vector<float> > norm_results(m_nPixels, std::vector<float>(m_nLabels));
      // tree_count is stored in node_in_trees[i].size()


      // for each node
      for (unsigned int i=0; i<node_in_subgraphs.size(); i++)
      {
        for (int t=0; t < node_in_subgraphs[i].size(); t++)
        {
          int subgraph_no = node_in_subgraphs[i][t][0];
          int subgraph_node_no = node_in_subgraphs[i][t][1];
          int label_no = mrf_vec[subgraph_no]->getLabel(subgraph_node_no);
          // no need to handle flipping as Flipper does it internally

          norm_results[i][label_no]++;
        }
      }

      // for each node
      for (unsigned int i=0; i<node_in_subgraphs.size(); i++)
      {
        for (unsigned int d=0; d<m_nLabels; d++)
        {
          norm_results[i][d] /= node_in_subgraphs[i].size();
        }
      }

      e_subgraphs_total_old = e_subgraphs_total;



      std::vector<int> mismatched_nodes;
      int num_mismatch = 0;

      for (int ii=0; ii<m_nPixels; ii++)
      {
        setLabel(ii, -1); // Label allows -ve numbers.
      }
      // for each node


      for (unsigned int ii=0; ii<node_in_subgraphs.size(); ii++)
      {
        int current_label = getLabel(ii);
        for (int t=0; t < node_in_subgraphs[ii].size(); t++)
        {
          int subgraph_no = node_in_subgraphs[ii][t][0];
          int subgraph_node_no = node_in_subgraphs[ii][t][1];
          int label_no = mrf_vec[subgraph_no]->getLabel(subgraph_node_no);
          // no need to handle flipping

          if (current_label == -1)
            current_label = label_no;
          else if (current_label != label_no)
          {
            changed = 1;
            num_mismatch++;
            mismatched_nodes.push_back(ii);
            break; // one mismatch per node.
            // printf(" Mismatched labels Graph node %d Tree No %d Tree node %d\n", i, tree_no, tree_node_no);
          }
        }
        if (current_label == -1)
          printf(" Serious representation error \n");
      }



      // std::vector<double> temp_primal_energies(trees_edge_list.size(), 0.0);
      double temp_best_primal = 0.0;
      double temp_best_primal_index = -1;
      // here we are assuming all trees are flippers generated from spanning trees.
      // for all subgraphs
      for (int h=0; h<subgraphs_edge_list.size(); h++)
      {
        // printf(" Tree %d \n", h);
        // we assume all nodes of graph are in flipper graph (created from spanning trees)
        for (int s=0; s<m_nPixels; s++)
        {
          int label_no = mrf_vec[h]->getLabel(s);
          int gnode = subgraph_node_to_graph_node[h][s];
          setLabel(gnode, label_no);
          // printf(" node %d label %d \n", gnode, label_no);
        }

        double primalv = totalEnergy();

        if (h == 0)
        {
          temp_best_primal = primalv;
          temp_best_primal_index = 0;
        } else {

          if (temp_best_primal > primalv)
          {
            temp_best_primal = primalv;
            temp_best_primal_index = h;
          }
        }
      }

      if (temp_best_primal_index == -1)
      {
        printf(" primal selecting error \n");
        exit(1);
      }


      for (int s=0; s<m_nPixels; s++)
      {
        int label_no = mrf_vec[temp_best_primal_index]->getLabel(s);
        int gnode = subgraph_node_to_graph_node[temp_best_primal_index][s];
        setLabel(gnode, label_no);
      }




      // printf("**** Finished Setting the labels ****\n");
      float primal_energy = temp_best_primal; // totalEnergy(); // this will use fnCost directly
      if (iter_no == 1)
      {
        best_primal_energy = primal_energy;
        for (unsigned int iii=0; iii<m_nPixels; iii++)
        {
          glabels[iii] = getLabel(iii);
        }
        gmismatchno = num_mismatch;
      }
      else
      {
        if (primal_energy < best_primal_energy)
        {
          best_primal_energy = primal_energy;
          for (unsigned int iii=0; iii<m_nPixels; iii++)
          {
            glabels[iii] = getLabel(iii);
          }
          gmismatchno = num_mismatch;
        }

        if (primal_energy == best_primal_energy && gmismatchno > num_mismatch)
        {
          best_primal_energy = primal_energy;
          for (unsigned int iii=0; iii<m_nPixels; iii++)
          {
            glabels[iii] = getLabel(iii);
          }
          gmismatchno = num_mismatch;
        }


      }

      clock_t finishss = clock();
      float timess = (float) (((double)(finishss - start)) / CLOCKS_PER_SEC);

      timess += timeoffset;

      // clock_t overflows after 2000 seconds and so need to restart it and store value
      // temporarily in timeoffset
      if (finishss < 0 && start >= 0)
      {
        timess = (float) (((double)(LONG_MAX + finishss) + double(LONG_MAX - start)) / (double)CLOCKS_PER_SEC);
        timeoffset += timess;
        timess = timeoffset;
        start = clock();
      }
      if (finishss > 0 && start < 0)
      {
        timess = (float) (((double)(finishss) + double(- start)) / (double)CLOCKS_PER_SEC);
        timeoffset += timess;
        timess = timeoffset;
        start = clock();
      }



      //printf("**PE = %f****DE = %f**\n", primal_energy, e_trees_total);
      // if (iter_no % 20 == 0)
        printf(" %d : mismatch = %d **  PE = %f  DE = %f ( %f secs)", iter_no, num_mismatch, primal_energy, e_subgraphs_total, timess);
      // if (iter_no % 20 == 0)
        printf("\n");

      iter_no++;
      if (changed == 0 || iter_no > max_iter || fabs(primal_energy - e_subgraphs_total) < getToleranceParam())
        break;

      // computing the norm of the gradient
      float norm_grad = 0;

      for (unsigned int ii=0; ii<mismatched_nodes.size(); ii++)
      {
         int i = mismatched_nodes[ii];

         for (unsigned int t=0; t < node_in_subgraphs[i].size(); t++)
         {
           int subgraph_no = node_in_subgraphs[i][t][0];
           int subgraph_node_no = node_in_subgraphs[i][t][1];
           int label_no = mrf_vec[subgraph_no]->getLabel(subgraph_node_no);
           // no need to handle flipping

           for (unsigned int d=0; d<m_nLabels; d++)
           {
             float gradient = -norm_results[i][d];
             if (d == label_no)
               gradient += 1;
             // else x_p_d is 0
             norm_grad += gradient*gradient;
           }
         }
       }



      // for each node of each tree, change datacost term
      // as per gradient.

      // printf("**** Adding gradient to datacost terms for next iteration ****\n");

      std::vector<int> subgraphs_needed_to_be_updated(num_subgraphs, 0);


      for (unsigned int ii=0; ii<mismatched_nodes.size(); ii++)
      {
        int i = mismatched_nodes[ii];
        // printf("**** i = %d \n", i);

        for (unsigned int t=0; t < node_in_subgraphs[i].size(); t++)
        {
          int subgraph_no = node_in_subgraphs[i][t][0];
          subgraphs_needed_to_be_updated[subgraph_no] = 1;
          int subgraph_node_no = node_in_subgraphs[i][t][1];
          int label_no = mrf_vec[subgraph_no]->getLabel(subgraph_node_no);


          // This x_i_treeno(label_no) = 1 all else x_i_treeno are 0
          // norm_results[i][label_no]++;
          // MaxProdBPTreeDD* mrf_temp = (MaxProdBPTreeDD*)mrf_vec[tree_no];
          int idx;

          for (unsigned int d=0; d<m_nLabels; d++) // d is an iterating index here
          {
            int di = d; // di is used to get the datacost term of the graph

            idx = subgraph_node_no*m_nLabels + di;
            float theta_p_d = mrf_vec[subgraph_no]->getEnergyFunction()->m_dataCost->getDataCostRaw(idx);

            float gradient = -norm_results[i][d];

            if (d == label_no)
              gradient += 1;
            // else x_p_d is 0

            // 0.00000000001; //
            gradient *= gamma*(best_primal_energy - e_subgraphs_total)/norm_grad;

            ((FlipperForDD*)mrf_vec[subgraph_no])->updateFlipperDataStructures(idx, subgraph_node_no, di, theta_p_d + gradient);

          }
        }
      }

      // alpha = 10.0/iter_no;
      // alpha *= 0.9;
      // alpha = 0.1/sqrt(double(iter_no));
    }



    //for (int i=0; i<m_nPixels; i++)
    //{
    //  setLabel(i, -1); // Label allows -ve numbers.
    //}

    // using best mismatch result
    for (unsigned int i=0; i<m_nPixels; i++)
    {
      setLabel(i, glabels[i]);
    }

    printf(" Mismatches : %d \n", gmismatchno);


    for (unsigned int t=0; t<num_subgraphs; t++)
    {
      MRF::CostVal* ptr = D_a[t];
      delete[] ptr;
    }

  }

}







