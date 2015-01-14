/*
 * DualDecompositionSubGraph.cpp
 *
 *  Created on: Dec 28, 2011
 *      Author: bhole
 */



#include <time.h>
#include <limits.h>
#include <math.h>
#include <queue>

#include "DualDecompositionSubGraph.h"
#include "GCoptimization.h"



DualDecompositionSubGraph::DualDecompositionSubGraph(int width, int height, int nLabels,EnergyFunction *eng):DualDecomposition(width,height,nLabels,eng)
{
}
DualDecompositionSubGraph::DualDecompositionSubGraph(int nPixels, int nLabels,EnergyFunction *eng):DualDecomposition(nPixels,nLabels,eng)
{
}

DualDecompositionSubGraph::~DualDecompositionSubGraph()
{
}

void DualDecompositionSubGraph::setSubmodularFunction(MRF::SmoothCostGeneralFn subfnCost)
{
  m_ssmoothFn = subfnCost;
}

void DualDecompositionSubGraph::setNonSubmodularFunction(MRF::SmoothCostGeneralFn nonsubfnCost)
{
  m_nssmoothFn = nonsubfnCost;
}

void DualDecompositionSubGraph::setflippedinfo(std::vector<int> **flippedinfo)
{
  flippedinfo_ptop = flippedinfo;
}

void DualDecompositionSubGraph::setNodeMappingStructure(
    std::vector<int> **sg_node_to_g_node)
{
  subgraph_node_to_graph_node_ptop = sg_node_to_g_node;
}

/*

#include "MaxProdBPTreeDD.h"
// #include "regions-new.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN






void DualDecomposition::extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
                                        std::vector<std::vector<int> > *edge_list)
{
  int edgenumber = 0;

  std::vector<int> temp_list(2,-1);

  for (int i = 0; i < m_nPixels; i++)
  {
    for (int p = 0; p < m_nNeighbors[i]; p++)
    {
      int q=findNeighborWithIdx(i, p);
      if (q >= i)
      {
        // store [other node, edge number]
        temp_list[0] = q;
        temp_list[1] = edgenumber;
        (*node_edge_list)[i].push_back(temp_list);

        // store [other node, edge number]
        temp_list[0] = i;
        (*node_edge_list)[q].push_back(temp_list);

        // store [first node, second node]
        temp_list[1] = q;
        edge_list->push_back(temp_list);

        edgenumber++;
      } // else it means it was put in already
    }
  }
}

int DualDecomposition::getEdgeForNodes(const int &start_node,  const int &neighbour,
                                        const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list)
{
  for (unsigned int i=0; i<graph_node_edge_list[start_node].size(); i++)
  {
    if (graph_node_edge_list[start_node][i][0] == neighbour)
    {
      return graph_node_edge_list[start_node][i][1];
    }
  }
  return -1; //fail
}

// 0 not prefered
// 1 prefered
void DualDecomposition::extractPrimTreeEdges(
    const int &prefer_not_used_edges,
    const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *global_used_edges,
    std::vector<std::vector<int> > *trees_edge_list)
{

  std::vector<int> tree_edge_list(m_nPixels-1, -1);
  // not seen nodes are 0, seen are 1
  std::vector<int> seen_nodes(m_nPixels, 0);
  int num_edges = 0;

  int start_node = -1;

  if (prefer_not_used_edges==1)
  {
    for (unsigned int e=0; e < global_used_edges->size(); e++)
    {
      if ((*global_used_edges)[e] == 0)
      {
        start_node = graph_edge_list[e][0];
        break;
      }
    }
  }
  // Even for a preferred case, if all edges are used, start node can be -1
  if (start_node == -1)
  {
    // randomly pick a node to start working on
    start_node = floor((double(rand())/ RAND_MAX)*(m_nPixels));
    if (start_node>=m_nPixels)
      start_node = m_nPixels - 1;
  }


  seen_nodes[start_node] = 1;
  //printf(" graph_edge_list.size() :%d ", graph_edge_list.size());
  // printf(" start node : %d \n", start_node);
  //printf(" graph_edge_list.size() :%d ", graph_edge_list.size());


  // no edge is currently used for the tree
  std::vector<int> tree_used_edges(graph_edge_list.size(), 0);

  std::vector<int> candidate_list;
  // put all its neighbors in candidate list

  //printf(" neighbor :%d \n", m_nNeighbors[start_node]);

  for (int p = 0; p < m_nNeighbors[start_node]; p++)
  {
    int q = findNeighborWithIdx(start_node, p);
    int edge = getEdgeForNodes(start_node,  q, graph_node_edge_list);
    //printf(" edge :%d ", edge);
    if (edge == -1)
    {
       printf(" Edge not found \n");
       exit(1);
    }
    candidate_list.push_back(edge);
  }

  //printf(" candidate_list :%d ", candidate_list.size());
  assert(candidate_list.size()>0);

  while(num_edges < m_nPixels-1)
  {
    int candidate_edge = -1;
    int candidate_edge_idx = -1;

    if (prefer_not_used_edges==1)
    {
      for (unsigned int i=0; i<candidate_list.size(); i++)
      {
        if ((*global_used_edges)[candidate_list[i]]==0)
        {
          candidate_edge = candidate_list[i];
          candidate_edge_idx = i;
          break;
        }
      }
    }
    // if all edges seem to have been used if we used prefered_not_seen_edges
    // then we randomly just pick one edge.
    if (candidate_edge == -1)
    {
      int idx_c = floor((double(rand())/ RAND_MAX)*(candidate_list.size()));
      if (idx_c >= candidate_list.size())
        idx_c = candidate_list.size() - 1;
      assert(idx_c>=0);
      candidate_edge = candidate_list[idx_c];
      candidate_edge_idx = idx_c;
    }

    // printf("\n candidate_edge : %d ", candidate_edge);
    // printf(" candidate_edge_idx : %d ", candidate_edge_idx);

    int current = -1;
    int other = -1;
    if (seen_nodes[graph_edge_list[candidate_edge][0]] == 1)
    {
      current = graph_edge_list[candidate_edge][0];
      other = graph_edge_list[candidate_edge][1];
    } else {
      current = graph_edge_list[candidate_edge][1];
      other = graph_edge_list[candidate_edge][0];
    }

    if (seen_nodes[other] == 1)
    {
      // then you can't use this edge because it will form a cycle.
      // Note that one of the seen_nodes will already be 1;
    } else {
      // use the edge
      // printf(" Node: %d ", other);
      seen_nodes[other] = 1;
      tree_used_edges[candidate_edge] = 1;
      (*global_used_edges)[candidate_edge] = 1;
      tree_edge_list[num_edges] = candidate_edge;
      num_edges++;

      // add neighbors of other in the candidate list
      for (int p = 0; p < m_nNeighbors[other]; p++)
      {
        int q = findNeighborWithIdx(other, p);
        if (seen_nodes[q] != 1)
        {
          int edge = getEdgeForNodes(other,  q, graph_node_edge_list);
          if (edge == -1)
          {
             printf(" Edge not found \n");
             exit(1);
          }
          candidate_list.push_back(edge);
        }
      }
    }

    // remove this edge from candidatelist
    for (unsigned int i=candidate_edge_idx+1; i < candidate_list.size(); i++)
    {
      candidate_list[i-1] = candidate_list[i];
    }
    int size = candidate_list.size();
    candidate_list.resize(size-1);
  }

  trees_edge_list->push_back(tree_edge_list);

}


void DualDecomposition::extractDFSTreeMaxLimitNonRepeatableEdges(
    const int &max_edge_limit,
    const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *global_used_edges,
    std::vector<std::vector<int> > *trees_edge_list)
{

  int max_edge_limit_curr = max_edge_limit;
  if (max_edge_limit_curr > m_nPixels-1)
  {
    max_edge_limit_curr = m_nPixels - 1;
  }

  std::vector<int> tree_edge_list;

  // not seen nodes are 0, seen are 1
  std::vector<int> seen_nodes(m_nPixels, 0);
  int num_edges = 0;

  int start_node = -1;

  for (unsigned int e=0; e < global_used_edges->size(); e++)
  {
    if ((*global_used_edges)[e] == 0)
    {
      start_node = graph_edge_list[e][0];
      break;
    }
  }

  if (start_node == -1)
  {
    // all edges are covered.
    return;
  }

  seen_nodes[start_node] = 1;

  // no edge is currently used for the tree
  std::vector<int> tree_used_edges(graph_edge_list.size(), 0);

  std::vector<int> candidate_list;
  // put all its neighbors in candidate list

  for (int p = 0; p < m_nNeighbors[start_node]; p++)
  {
    int q = findNeighborWithIdx(start_node, p);
    int edge = getEdgeForNodes(start_node,  q, graph_node_edge_list);
    if (edge == -1)
    {
       printf(" Edge not found \n");
       exit(1);
    }
    if ((*global_used_edges)[edge]==0)
    {
      candidate_list.push_back(edge);
    }
  }

  assert(candidate_list.size()>0);

  int current = start_node;
  int other = -1;

  while(num_edges < max_edge_limit_curr)
  {
    int candidate_edge = -1;
    int candidate_edge_idx = -1;

    for (int i=candidate_list.size()-1; i>=0; i--)
    {
      if ((*global_used_edges)[candidate_list[i]]==0)
      {
        candidate_edge = candidate_list[i];
        candidate_edge_idx = i;
        break;
      }
    }

    if (candidate_edge == -1)
    {
      // all edges are covered or used
      break;
    }

    // printf(" vv : %d ", candidate_list.size());

    if (other == -1 && current == -1)
    {
      if (seen_nodes[graph_edge_list[candidate_edge][0]] == 1)
      {
        current = graph_edge_list[candidate_edge][0];
        other = graph_edge_list[candidate_edge][1];
      } else {
        current = graph_edge_list[candidate_edge][1];
        other = graph_edge_list[candidate_edge][0];
      }
    } else {
      if (other == -1 && current != -1)
      {
        int curr_t = graph_edge_list[candidate_edge][0];
        int other_t = graph_edge_list[candidate_edge][1];

        if (curr_t == current)
          other = other_t;
        else if (other_t == current)
          other = curr_t;
        else {
          // the edge belongs to a parent or grand parent of previous node
          if (seen_nodes[curr_t] == 1)
          {
            current = curr_t;
            other = other_t;
          } else {
            current = other_t;
            other = curr_t;
          }
        }
      }
    }

    if (seen_nodes[other] == 1)
    {
      // then you can't use this edge because it will form a cycle.
      // Note that one of the seen_nodes will already be 1;
      current = -1;
      other = -1;
    } else {
      // use the edge
      // printf(" Node: %d ", other);
      seen_nodes[other] = 1;
      tree_used_edges[candidate_edge] = 1;
      (*global_used_edges)[candidate_edge] = 1;
      tree_edge_list.push_back(candidate_edge);
      num_edges++;

      // add neighbors of other in the candidate list
      for (int p = 0; p < m_nNeighbors[other]; p++)
      {
        int q = findNeighborWithIdx(other, p);
        if (seen_nodes[q] != 1)
        {
          int edge = getEdgeForNodes(other,  q, graph_node_edge_list);
          if (edge == -1)
          {
             printf(" Edge not found \n");
             exit(1);
          }
          if ((*global_used_edges)[edge]==0)
          {
            candidate_list.push_back(edge);
          }
        }
      }

      current = other;
      other = -1;
    }

    // remove this edge from candidatelist
    for (unsigned int i=candidate_edge_idx+1; i < candidate_list.size(); i++)
    {
      candidate_list[i-1] = candidate_list[i];
    }
    int size = candidate_list.size();
    candidate_list.resize(size-1);
  }

  if (tree_edge_list.size() != 0)
    trees_edge_list->push_back(tree_edge_list);

}





int areThereUnusedEdges(const std::vector<int> &global_used_edges)
{
  int unused = 0;
  for (unsigned int i=0; i < global_used_edges.size(); i++)
  {
    if (global_used_edges[i] == 0)
    {
      unused = 1;
      break;
    }
  }
  return unused;
}


*/


void DualDecompositionSubGraph::getSubGraphNodeNo(const std::vector<std::vector<int> > &node_in_subgraphs_g,
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



void DualDecompositionSubGraph::getWhatSubGraphsNodeIsUsedIn(const std::vector<std::vector<int> > &subgraphs_edge_list,
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





void DualDecompositionSubGraph::dumpSubGraphs(const std::vector<std::vector<int> > &subgraphs_edge_list,
               const std::vector<std::vector<int> > &edge_list)
{
  for (unsigned int i=0; i < subgraphs_edge_list.size(); i++)
  {
    printf("\n SubGraph : %d SubGraph size : %d\n", i, subgraphs_edge_list[i].size());
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



// 1 for submodular, 0 for non submodular
int DualDecompositionSubGraph::checkSubmodularityEdge(const int &node1, const int &node2)
{
  if ( m_grid_graph )
  {
    printf(" Grid graphs don't work \n");
    exit(1);
  }

  if (! m_grid_graph) // general graph
  {
    if ( m_smoothType != FUNCTION  )
    {
      printf(" Fixed array shouldn't work for now \n");
      exit(1);
    }
    else // smooth is a function
    {
      MRF::CostVal E00 = m_smoothFn(node1, node2, 0, 0);
      MRF::CostVal E01 = m_smoothFn(node1, node2, 0, 1);
      MRF::CostVal E10 = m_smoothFn(node1, node2, 1, 0);
      MRF::CostVal E11 = m_smoothFn(node1, node2, 1, 1);

      if (E00 + E11 <= E10 + E01)
        return 1;
      else
        return 0;
    }
  }
}

// need to make sure flipped subgraph does not have a graph that will generate the flipping inconsistently
// ie a flipped node cannot be connected to two unflipped nodes that are also neighboring each other

void DualDecompositionSubGraph::extractSubGraphEdges(
    const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *global_used_edges,
    std::vector<std::vector<int> > *subgraphs_edge_list,
    std::vector<int> *subgraph_type)
{

  std::vector<int> submod_edge_list;
  std::vector<int> nonsubmod_edge_list;

  //std::vector<int> flipped_nodes;

  for (unsigned int e=0; e < global_used_edges->size(); e++)
  {
    if ((*global_used_edges)[e] == 0)
    {
      int first_node = graph_edge_list[e][0];
      int second_node = graph_edge_list[e][1];

      if (checkSubmodularityEdge(first_node, second_node) == 1) // submodular
      {
        submod_edge_list.push_back(e);
      } else
      { // non submodular
        nonsubmod_edge_list.push_back(e);
      }

      (*global_used_edges)[e] = 1;

    }
  }

  if (submod_edge_list.size() > 0)
  {
    subgraphs_edge_list->push_back(submod_edge_list);
    subgraph_type->push_back(0); //  graph is submodular
  }
  if (nonsubmod_edge_list.size() > 0)
  {
    subgraphs_edge_list->push_back(nonsubmod_edge_list);
    subgraph_type->push_back(1); //  graph is non-submodular
  }

}



void populateNodeNeighborList(const std::vector<int> &subgraph_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<std::vector <int> > *node_neighbor)
{
  for (unsigned int e = 0; e < subgraph_edge_list.size(); e++)
  {
    int edge = subgraph_edge_list[e];

    int node1 = edge_list[edge][0];
    int node2 = edge_list[edge][1];

    int already_in = 0;
    for (int i=0; i < (*node_neighbor)[node1].size(); i++)
    {
      if ((*node_neighbor)[node1][i] == node2)
      {
        already_in = 1;
        break;
      }
    }
    if (already_in == 0)
      (*node_neighbor)[node1].push_back(node2);

    already_in = 0;
    for (int i=0; i < (*node_neighbor)[node2].size(); i++)
    {
      if ((*node_neighbor)[node2][i] == node1)
      {
        already_in = 1;
        break;
      }
    }
    if (already_in == 0)
      (*node_neighbor)[node2].push_back(node1);
  }
}



void printNodeNeighborList(const std::vector<std::vector <int> > &node_neighbor)
{
  for (unsigned int i = 0; i < node_neighbor.size(); i++)
  {
    printf(" Node : %d  Neighbors :: ", i );
    for (unsigned int j = 0; j < node_neighbor[i].size(); j++)
    {
      printf(" %d ", node_neighbor[i][j]);
    }
    printf("\n");
  }
}


// mark them in node_in_subgraphs
void DualDecompositionSubGraph::markFlippedNodes(
    const std::vector<std::vector<int> > &subgraphs_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &subgraph_type,
    std::vector<std::vector<std::vector<int> > >*node_in_subgraphs)
{
  for (unsigned int t = 0; t < subgraphs_edge_list.size(); t++)
  {
    if (subgraph_type[t] == 1) // is non-submodular
    {
      // traverse this subgraph t

      // creating a node neighbor list
      std::vector<std::vector <int> > node_neighbor(m_nPixels, (std::vector <int>)0);
      populateNodeNeighborList(subgraphs_edge_list[t], edge_list, &node_neighbor);

      // printNodeNeighborList(node_neighbor);

      std::queue<int> flippedlist;
      std::vector<int> checked(m_nPixels, 0);
      // because there can be forests, we need to go thru all nodes

      for (int i=0; i<m_nPixels; i++)
      {
        if (checked[i] == 0 && node_neighbor[i].size() > 0) // this allows start of a new forest
        {
          checked[i] = 1;
          for (int nn=0; nn < node_neighbor[i].size(); nn++)
          {
            int neighbor = node_neighbor[i][nn];
            if (checked[neighbor] == 0)
            {
              flippedlist.push(neighbor);
              checked[neighbor] = 1;
            }
          }

          while (!flippedlist.empty())
          {
            int nodem = flippedlist.front();
            flippedlist.pop();
            //printf(" %d ", nodem);

            // put this node as flipped in node_in_subgraphs
            for (int kk=0; kk<(*node_in_subgraphs)[nodem].size(); kk++)
              if ((*node_in_subgraphs)[nodem][kk][0] == t)
                (*node_in_subgraphs)[nodem][kk][2] = 1;

            // these are the non-flipped nodes
            for (int ns=0; ns < node_neighbor[nodem].size(); ns++)
            {
              int noden = node_neighbor[nodem][ns];
              // printf(" *%d* ", noden);

              if (checked[noden] == 0)
              {
                checked[noden] = 1;
                // put neighbors of this node in list
                for (int nt=0; nt < node_neighbor[noden].size(); nt++)
                {
                  int neighbor = node_neighbor[noden][nt];
                  if (checked[neighbor] == 0)
                  {
                    flippedlist.push(neighbor);
                    checked[neighbor] = 1;
                  }
                }
              }

            }
          }

        }
      }
    }
  }
}


// the node number order is maintained.
void DualDecompositionSubGraph::computeDataCostForSubGraphs(
    const std::vector<std::vector<int> > &subgraphs_edge_list,
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &subgraph_type,
    std::vector<std::vector<std::vector<int> > > *node_in_subgraphs,
    std::vector<MRF::CostVal*> *D_a,
    std::vector<std::vector<int> > *subgraph_node_to_graph_node,
    std::vector<DataCost*> *dcost_vec,
    std::vector<std::vector<int> > *flipstructure)
{

  assert(m_nLabels == 2);

  std::vector<int> dsiIndex(subgraphs_edge_list.size(), 0);

  // node_in_subgraphs gives all the subgraphs that a node belongs to
  for (unsigned int i=0; i < node_in_subgraphs->size(); i++)
  {
    assert((*node_in_subgraphs)[i].size() != 0);
    for (unsigned int t = 0; t < (*node_in_subgraphs)[i].size(); t++)
    {
      int subgraph_num = (*node_in_subgraphs)[i][t][0];
      (*subgraph_node_to_graph_node)[subgraph_num].push_back(i);
      (*flipstructure)[subgraph_num].push_back(0);
      int subgraph_node_no = (*subgraph_node_to_graph_node)[subgraph_num].size() - 1;
      (*node_in_subgraphs)[i][t][1] = subgraph_node_no;
      if (subgraph_type[subgraph_num] == 0 || (*node_in_subgraphs)[i][t][2] == 0) // submodular or non-flip node
      {
        for (int d = 0; d < m_nLabels; d++)
        {
          MRF::CostVal* ptr = (*D_a)[subgraph_num] + dsiIndex[subgraph_num];
          *ptr = nodeArray[i].localEv[d] / (*node_in_subgraphs)[i].size();
          dsiIndex[subgraph_num]++;
        }
      } else // non-submodular and flip node
      { // datacost term values are reversed
        for (int d = m_nLabels-1; d >= 0; d--)
        {
          MRF::CostVal* ptr = (*D_a)[subgraph_num] + dsiIndex[subgraph_num];
          *ptr = nodeArray[i].localEv[d] / (*node_in_subgraphs)[i].size();
          dsiIndex[subgraph_num]++;
        }
        // update flip structure because this is flip node
        (*flipstructure)[subgraph_num][subgraph_node_no] = 1;
      }
    }
  }

  for (unsigned int t = 0; t < subgraphs_edge_list.size(); t++)
  {
    (*dcost_vec)[t] = new DataCost ((*D_a)[t]);
  }
}




void DualDecompositionSubGraph::optimizeAlg(int nIterations)
{
  if (m_grid_graph)
  {
    printf(" Dual Decomposition Tree does not handle grid graphs \n");
    exit(1);
  } else {
    // for non-grid graphs

    // Until number of iterations or convergence.
    // 1. use current graph to create subgraphs
    // 2. Split datacost and smoothcost terms for subgraphs
    // 3. use Graph cuts for each subgraph (after flipping trick) to do inference
    // 4. update parameters using projected sub-gradients.

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

    // start timer to compute how much time it takes to create subgraphs
    // and copy datacost, smoothcost etc
    //clock_t startss = clock();



    // [subgraph]
    std::vector<int> subgraph_type; // 0 for submodular and 1 for non-submodular
    extractSubGraphEdges(node_edge_list, edge_list, &global_used_edges, &subgraphs_edge_list, &subgraph_type);
    // you will need to make sure that the non-submodular sub graph does not contain cycles of odd nodes or edges
    // else it means flipping trick won't work

    // after calling all subgraphs dump edges for each subgraph for debugging purpose.
    // dumpSubGraphs(subgraphs_edge_list, edge_list);

    const int num_subgraphs = subgraphs_edge_list.size();

    // It is important to note the node numbering used for subgraphs.
    // ie what node numbers in the subgraph correspond to what node numbers
    // of the full graph
    // same for the edges.
    // hence we have mapping functions for both these.
    // the node number mapping is got from computeDataCostForTree

    // flip variable is used for constructing the subgraph that does not follow submodularity
    // this data structure is [graph_node_no][subgraph][0:subgraph_no, 1:subgraph_node_no 2:flipped]
    std::vector<std::vector<std::vector<int> > > node_in_subgraphs(m_nPixels, (std::vector<std::vector<int> >)0);
    // this structure counts number of nodes in each subgraph
    std::vector<int> num_nodes_in_subgraphs(subgraphs_edge_list.size(), 0);
    getWhatSubGraphsNodeIsUsedIn(subgraphs_edge_list, edge_list, &node_in_subgraphs, &num_nodes_in_subgraphs);

    // printf("%d %d\n", num_nodes_in_subgraphs[0], num_nodes_in_subgraphs[1]);


    markFlippedNodes(subgraphs_edge_list, edge_list, subgraph_type, &node_in_subgraphs);

    for (unsigned int i=0; i < node_in_subgraphs.size(); i++)
    {
      // make sure that each node is used at least in one tree.
      assert(node_in_subgraphs[i].size() != 0);
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
    std::vector<std::vector<int> > flipstructure(num_subgraphs, (std::vector<int>)0);

    printf(" Computing DataCost for trees ... \n");
    computeDataCostForSubGraphs(subgraphs_edge_list, edge_list, subgraph_type, &node_in_subgraphs, &D_a, &subgraph_node_to_graph_node, &dcost_vec, &flipstructure);
    // datacost with flipped nodes have datacost terms flipped.
    // so we need to remember that flipping while changing the parameters

    SmoothnessCost* sscost_vec = new SmoothnessCost((MRF::SmoothCostGeneralFn) m_ssmoothFn); // submodular
    SmoothnessCost* nsscost_vec = new SmoothnessCost((MRF::SmoothCostGeneralFn) m_nssmoothFn); // nonsubmodular


    // create mrfs of Expansion for each subgraph we have.
    for (int t=0; t < num_subgraphs; t++)
    {
      if (subgraph_type[t] == 0)
      {
        energy_vec[t] = new EnergyFunction(dcost_vec[t], sscost_vec);
      } else
      {
        energy_vec[t] = new EnergyFunction(dcost_vec[t], nsscost_vec);
      }

      mrf_vec[t] = new Expansion(subgraph_node_to_graph_node[t].size(), m_nLabels, energy_vec[t]);

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
    // 200
    int max_iter = 10; // refers to infer from subgraph not gradient descent times
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
        *flippedinfo_ptop = &flipstructure[t];
        *subgraph_node_to_graph_node_ptop = &subgraph_node_to_graph_node[t];
        E = mrf_vec[t]->totalEnergy();
        tot_time_t= 0;
        mrf_vec[t]->optimize(1, time_t);
        E = mrf_vec[t]->totalEnergy();
        tot_time_t = tot_time_t + time_t ;

/*        for (int ix=0; ix<num_nodes_in_subgraphs[t]; ix++)
        {
          unsigned char rowvalue = mrf_vec[t]->getLabel(ix);
          printf("Tree : %d  Node : %d  state : %d \n", t, ix, (int)rowvalue);
        }
*/

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
          // handle flipping
          if (subgraph_type[subgraph_no] == 1 && flipstructure[subgraph_no][subgraph_node_no] == 1)
            label_no = 1 - label_no;

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
//      if (iter_no == 1)
//      {
//        changed = 1;
//      } else
      {
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
            // handle flipping
            if (subgraph_type[subgraph_no] == 1 && flipstructure[subgraph_no][subgraph_node_no] == 1)
              label_no = 1 - label_no;
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
          else
            setLabel(ii, current_label);
        }
      }

      // save labels of best mismatch solution
      // since iteration 1 starts with all tree nodes assigned to same value 0
      // and then calling mrf_vec[t]->optimize(1, time_t);, mismatch should not really be 0
      // unless the even breaking of parameter values was good.
      // if (num_mismatch < gmismatchno) // capturing the best solution in glabels
      // {
      //    for (unsigned int iii=0; iii<m_nPixels; iii++)
      //    {
      //      glabels[iii] = getLabel(iii);
      //    }
      //    gmismatchno = num_mismatch;
      //  }


      // printf("**** Finished Setting the labels ****\n");
      float primal_energy = totalEnergy();
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
           // handle flipping
           if (subgraph_type[subgraph_no] == 1 && flipstructure[subgraph_no][subgraph_node_no] == 1)
             label_no = 1 - label_no;
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

          if (subgraph_type[subgraph_no] == 1 && flipstructure[subgraph_no][subgraph_node_no] == 1)
            label_no = 1 - label_no;

          // This x_i_treeno(label_no) = 1 all else x_i_treeno are 0
          // norm_results[i][label_no]++;
          // MaxProdBPTreeDD* mrf_temp = (MaxProdBPTreeDD*)mrf_vec[tree_no];
          int idx;

          for (unsigned int d=0; d<m_nLabels; d++) // d is an iterating index here
          {
            int di = d; // di is used to get the datacost term of the graph
            if (subgraph_type[subgraph_no] == 1 && flipstructure[subgraph_no][subgraph_node_no] == 1)
              di=1-di; // since the datacost terms are flipped for non-sub graph

            idx = subgraph_node_no*m_nLabels + di;
            float theta_p_d = mrf_vec[subgraph_no]->getEnergyFunction()->m_dataCost->getDataCostRaw(idx);

            float gradient = -norm_results[i][d];

            if (d == label_no)
              gradient += 1;
            // else x_p_d is 0

            // 0.00000000001; //
            gradient *= gamma*(best_primal_energy - e_subgraphs_total)/norm_grad;
            mrf_vec[subgraph_no]->getEnergyFunction()->m_dataCost->updateDataCostRaw(idx, theta_p_d + gradient);
          }
        }
      }

      // I don't have nodeArray in Expansion or GC and so don't need to update it.
      // but i need to double check that updateDataCostRaw fixes everything.
      /*
      for (int i=0; i < num_subgraphs; i++)
      {
        if (subgraphs_needed_to_be_updated[i] == 1)
        {
          //Expansion* mrf_temp = (Expansion*)mrf_vec[i];
          MaxProdBPTreeDD* mrf_temp = (MaxProdBPTreeDD*)mrf_vec[i];
          mrf_temp->updateData(mrf_vec[i]->getEnergyFunction()->m_dataCost->getDataCostArray());
        }
      }
      */

      // alpha = 10.0/iter_no;
      // alpha *= 0.9;
      // alpha = 0.1/sqrt(double(iter_no));
    }



    for (int i=0; i<m_nPixels; i++)
    {
      setLabel(i, -1); // Label allows -ve numbers.
    }

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




// for each node
// Use this only when we are sure trees have converged
/*
for (unsigned int i=0; i<node_in_trees.size(); i++)
{
  int current_label = getLabel(i);
  for (int t=0; t < node_in_trees[i].size(); t++)
  {
    int tree_no = node_in_trees[i][t][0];
    int tree_node_no = node_in_trees[i][t][1];
    int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);
    if (current_label == -1)
      current_label = label_no;
    else if (current_label != label_no)
    {
      printf(" Mismatched labels Graph node %d Tree No %d Tree node %d\n", i, tree_no, tree_node_no);
    }
  }
  if (current_label == -1)
    printf(" Serious representation error \n");
  else
    setLabel(i, current_label);
}
*/


