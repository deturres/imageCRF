/*
 * DualDecomposition.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: bhole
 */



#include <time.h>
#include <math.h>
#include <queue>

#include "DualDecomposition.h"
#include "MaxProdBPTreeDD.h"
// #include "regions-new.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN




DualDecomposition::DualDecomposition(int width, int height, int nLabels,EnergyFunction *eng):MRFWrapper(width,height,nLabels,eng)
{
  primaldualtolerance = 1e-10;
}
DualDecomposition::DualDecomposition(int nPixels, int nLabels,EnergyFunction *eng):MRFWrapper(nPixels,nLabels,eng)
{
  primaldualtolerance = 1e-10;
}

DualDecomposition::~DualDecomposition()
{

}


int DualDecomposition::areThereUnusedEdges(const std::vector<int> &global_used_edges)
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








int DualDecomposition::attemptAddEdge(const int curr_edge,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *nodes_in_what_tree,
    int *kruskal_forest_size,
    std::vector<std::vector <int> > *kruskal_forest,
    int *num_edges)
{

  int success = 1;

  int node1 = graph_edge_list[curr_edge][0];
  int node2 = graph_edge_list[curr_edge][1];

  if ((*nodes_in_what_tree)[node1] == -1 && (*nodes_in_what_tree)[node2] == -1)
  {
    // add new kruskal tree
    (*nodes_in_what_tree)[node1] = *kruskal_forest_size;
    (*nodes_in_what_tree)[node2] = *kruskal_forest_size;
    kruskal_forest->push_back((std::vector<int>)0);
    (*kruskal_forest)[*kruskal_forest_size].push_back(node1);
    (*kruskal_forest)[*kruskal_forest_size].push_back(node2);
    (*kruskal_forest_size)++;
    (*num_edges)++;
  } else if ((*nodes_in_what_tree)[node1] == -1 && (*nodes_in_what_tree)[node2] != -1)
  {
    // get tree number and add other node to that tree
    int tree_no = (*nodes_in_what_tree)[node2];
    (*nodes_in_what_tree)[node1] = tree_no;
    (*kruskal_forest)[tree_no].push_back(node1);
    (*num_edges)++;

  } else if ((*nodes_in_what_tree)[node1] != -1 && (*nodes_in_what_tree)[node2] == -1)
  {
    // get tree number and add other node to that tree
    int tree_no = (*nodes_in_what_tree)[node1];
    (*nodes_in_what_tree)[node2] = tree_no;
    (*kruskal_forest)[tree_no].push_back(node2);
    (*num_edges)++;

  } else if ((*nodes_in_what_tree)[node1] != -1 && (*nodes_in_what_tree)[node2] != -1
      && (*nodes_in_what_tree)[node1] != (*nodes_in_what_tree)[node2])
  {
    // merge the two trees
    int keep_tree_no = (*nodes_in_what_tree)[node1];
    int remove_tree_no = (*nodes_in_what_tree)[node2];

    for (int i=0; i<(*kruskal_forest)[remove_tree_no].size(); i++)
    {
      int node = (*kruskal_forest)[remove_tree_no][i];
      if (node == -1)
      {
        printf(" Some mistake in kruskal forest creation \n");
        exit(1);
      }
      (*nodes_in_what_tree)[node] = keep_tree_no;
      (*kruskal_forest)[remove_tree_no][i] = -1; // now that tree will have invalid nodes
    }
    (*num_edges)++;

  } else if ((*nodes_in_what_tree)[node1] != -1 && (*nodes_in_what_tree)[node2] != -1
      && (*nodes_in_what_tree)[node1] == (*nodes_in_what_tree)[node2])
  {
    // this will form a cycle in the tree so ignore this edge
    // or the edge is already included in the tree
    success = 0;
  }

  return success;

}


void printKruskalforest(const std::vector<std::vector <int> > &kruskal_forest,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> nodes_in_what_tree)
{
  for (int t=0; t<kruskal_forest.size(); t++)
  {
    for (int n=0; n<kruskal_forest[t].size(); n++)
    {
      int node = kruskal_forest[t][n];
      printf(" Tree : %d  Node %d \n", t, node );
    }
  }
}





void printTreeEdgeList(const std::vector<int> &tree_edge_list,
    const std::vector<std::vector<int> > &graph_edge_list)
{
  for (int e=0; e<tree_edge_list.size(); e++)
  {
    int node1 = graph_edge_list[tree_edge_list[e]][0];
    int node2 = graph_edge_list[tree_edge_list[e]][1];

    printf(" Edge %d  Node %d  Node %d\n", tree_edge_list[e], node1, node2 );
  }
}



void DualDecomposition::extractKruskalMaxLimitRepeatableEdges(
    const int &max_edge_limit,
    const std::vector<std::vector<std::vector<int> > > &graph_node_edge_list,
    const std::vector<std::vector<int> > &graph_edge_list,
    std::vector<int> *global_used_edges,
    std::vector<std::vector<int> > *trees_edge_list)
{

  const int global_num_nodes = graph_node_edge_list.size();


  int max_edge_limit_curr = max_edge_limit;
  if (max_edge_limit_curr > m_nPixels-1)
  {
    max_edge_limit_curr = m_nPixels - 1;
  }

  // number of edges in tree
  int num_edges = 0;
  // what edges does the final spanning tree have is stored in this structure
  std::vector<int> tree_edge_list;

  // create a priority edges list.
  // edges that are not already used have highest priority followed by all other edges in
  // the graph. // perhaps no need of storing other edges
  std::queue<int> priority_edges_list;

  for (unsigned int e=0; e < global_used_edges->size(); e++)
  {
    if ((*global_used_edges)[e] == 0)
    {
      priority_edges_list.push(e);
    }
  }

  // [tree_no [node_no]*]*
  std::vector<std::vector <int> > kruskal_forest;

  // this stores what kruskal tree a particular node belongs to currently
  // when two trees get connected, one tree needs to be destroyed after updating
  // the new kruskal_tree_no of the nodes merging.
  // [node_no kruskal_tree_no]*
  std::vector<int> nodes_in_what_tree(global_num_nodes, -1);


  int kruskal_forest_size = 0; // currently no tree

  // add as many edges from priority list as possible
  // condition is that if the edge has both the nodes from the same kruskal tree,
  // then it will form a cycle
  // otherwise you can add that edge so there are three cases
  // if none of the nodes is in any tree, then add new kruskal tree
  // else if one node is in one tree, other in no tree, add edge to that tree
  // else if one node is in one tree, other in another tree, you need to connect the two trees

  while (!priority_edges_list.empty())
  {
    int curr_edge = priority_edges_list.front();
    priority_edges_list.pop();

    // printf("Current edge priority %d\n", curr_edge);

    int success = attemptAddEdge(curr_edge, graph_edge_list, &nodes_in_what_tree, &kruskal_forest_size,
        &kruskal_forest, &num_edges);

    if (success == 1)
    {
      (*global_used_edges)[curr_edge] = 1;
      tree_edge_list.push_back(curr_edge);
    }
    // printf(" Printing kruskal forest \n");
    // printKruskalforest(kruskal_forest, graph_edge_list, nodes_in_what_tree);
    //printTreeEdgeList(tree_edge_list, graph_edge_list);

    if (num_edges >= max_edge_limit_curr)
      break;
  }



  // now just go thru all the edges in the global graph to make it a spanning tree

  for (int e=0; e<graph_edge_list.size(); e++)
  {
    if (num_edges >= max_edge_limit_curr)
      break;
    int success = attemptAddEdge(e, graph_edge_list, &nodes_in_what_tree, &kruskal_forest_size,
        &kruskal_forest, &num_edges);

    if (success == 1)
    {
      tree_edge_list.push_back(e);
      //printf("Current edge %d\n", e);
      //printTreeEdgeList(tree_edge_list, graph_edge_list);
    }

  }

  if (num_edges != max_edge_limit_curr)
  {
    printf(" Did not form a spanning tree, some problem in code \n");
    exit(1);
  }

  if (num_edges != tree_edge_list.size())
  {
    printf(" tree_edge_list and num_edges do not agree, some problem in code \n");
    exit(1);
  }

  // error checking step
  // check that all nodes (since it is spanning tree) have same tree number
  int main_tree_no = nodes_in_what_tree[0];

  for (int i=0; i<tree_edge_list.size(); i++)
  {
    int edge = tree_edge_list[i];

    int node1 = graph_edge_list[edge][0];
    int node2 = graph_edge_list[edge][1];

    if (nodes_in_what_tree[node1] != main_tree_no || nodes_in_what_tree[node2] != main_tree_no)
    {
      printf(" node kruskal tree consistency problem \n");
    }

  }


  if (tree_edge_list.size() != 0)
    trees_edge_list->push_back(tree_edge_list);
}

