/*
 * DualDecompositionTree.cpp
 *
 *  Created on: Dec 28, 2011
 *      Author: bhole
 */



#include <time.h>
#include <limits.h>
#include <math.h>

#include "DualDecompositionTree.h"
#include "MaxProdBPTreeDD.h"
// #include "regions-new.h"

#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN




DualDecompositionTree::DualDecompositionTree(int width, int height, int nLabels,EnergyFunction *eng):DualDecomposition(width,height,nLabels,eng)
{
}
DualDecompositionTree::DualDecompositionTree(int nPixels, int nLabels,EnergyFunction *eng):DualDecomposition(nPixels,nLabels,eng)
{
}

DualDecompositionTree::~DualDecompositionTree()
{
}


// 0 not prefered
// 1 prefered
void DualDecompositionTree::extractPrimTreeEdges(
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

  // no edge is currently used for the tree
  std::vector<int> tree_used_edges(graph_edge_list.size(), 0);

  std::vector<int> candidate_list;
  // put all its neighbors in candidate list

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






void dumpTrees(const std::vector<std::vector<int> > &trees_edge_list,
               const std::vector<std::vector<int> > &edge_list)
{
  for (unsigned int i=0; i < trees_edge_list.size(); i++)
  {
    printf("\n Tree : %d Tree size : %d\n", i, trees_edge_list[i].size());
    for (unsigned int j = 0; j < trees_edge_list[i].size(); j++)
    {
      int edge = trees_edge_list[i][j];
      if (edge == -1)
      {
        printf(" Error in edge ");
        exit(1);
      }
      printf("Edge: %d Node: %d Node: %d\n", edge, edge_list[edge][0], edge_list[edge][1]);
    }
  }
}



void getWhatTreesNodeIsUsedIn(const std::vector<std::vector<int> > &trees_edge_list,
		const std::vector<std::vector<int> > &edge_list,
		std::vector<std::vector<std::vector<int> > > *node_in_trees)
{
  const int num_trees = trees_edge_list.size();
  int edge_number;

  for (int t=0; t < num_trees; t++)
  {
    for (unsigned int e=0; e < trees_edge_list[t].size(); e++)
    {
      edge_number = trees_edge_list[t][e];
      int node1 = edge_list[edge_number][0];
      int node2 = edge_list[edge_number][1];

      // for node1
      int trees_size = (*node_in_trees)[node1].size();
      int found = 0;
      for (int f=0; f<trees_size; f++)
      {
        if (t == (*node_in_trees)[node1][f][0])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
      {
        std::vector<int> inn(2);
        inn[0] = t;
        inn[1] = -1; // tree node number will be filled in later function
        (*node_in_trees)[node1].push_back(inn);
      }

      // for node2
      trees_size = (*node_in_trees)[node2].size();
      found = 0;
      for (int f=0; f<trees_size; f++)
      {
        if (t == (*node_in_trees)[node2][f][0])
        {
          found = 1;
          break;
        }
      }
      if (found == 0)
      {
        std::vector<int> inn(2);
        inn[0] = t;
        inn[1] = -1; // tree node number will be filled in later function
        (*node_in_trees)[node2].push_back(inn);
      }
    }
  }

}

// the node number order is maintained.
void DualDecompositionTree::computeDataCostForTrees(
    const std::vector<std::vector<int> > &trees_edge_list,
		const std::vector<std::vector<int> > &edge_list,
		std::vector<std::vector<std::vector<int> > > &node_in_trees,
		std::vector<MRF::CostVal*> *D_a,
		std::vector<std::vector<int> > *tree_node_to_graph_node,
		std::vector<DataCost*> *dcost_vec)
{

  std::vector<int> dsiIndex(trees_edge_list.size(), 0);

  // node_in_trees gives all the trees that a node belongs to
  for (unsigned int i=0; i < node_in_trees.size(); i++)
  {
    assert(node_in_trees[i].size() != 0);
    for (unsigned int t = 0; t < node_in_trees[i].size(); t++)
    {
      int tree_num = node_in_trees[i][t][0];
      (*tree_node_to_graph_node)[tree_num].push_back(i);
      int tree_node_no = (*tree_node_to_graph_node)[tree_num].size() - 1;
      node_in_trees[i][t][1] = tree_node_no;
      for (int d = 0; d < m_nLabels; d++)
      {
        MRF::CostVal* ptr = (*D_a)[tree_num] + dsiIndex[tree_num];
        *ptr = nodeArray[i].localEv[d] / node_in_trees[i].size();
        dsiIndex[tree_num]++;
      }
    }
  }

  for (unsigned int t = 0; t < trees_edge_list.size(); t++)
  {
    (*dcost_vec)[t] = new DataCost ((*D_a)[t]);
  }
}


void getTreeNodeNo(const std::vector<std::vector<int> > &node_in_trees_g,
                   const int &tree_no, int *tree_node1)
{
  *tree_node1 = -1;
  for (unsigned int i = 0; i < node_in_trees_g.size(); i++)
  {
    if (tree_no == node_in_trees_g[i][0])
    {
      *tree_node1 = node_in_trees_g[i][1];
    }
  }
  assert(*tree_node1 != -1);
}







void DualDecompositionTree::optimizeAlg(int nIterations)
{
  if (m_grid_graph)
  {
    printf(" Dual Decomposition Tree does not handle grid graphs \n");
    exit(1);
  } else {
    // for non-grid graphs

    // Until number of iterations or convergence.
    // 1. use current graph to create multiple trees
    // 2. Split datacost and smoothcost terms for trees
    // 3. use MaxProdBPTree for each tree to do inference
    // 4. update parameters using projected sub-gradients.

    // If all edges while generating trees are not used, then additional trees will have to be made
    // Note that each tree will use all nodes.
    // Note if the graph is already a tree - you should really be using MaxProdBPTree

    // srand ( time(NULL) );
    //srand(0);

    // [node_no, [other_no, edge_no]*]*
    std::vector<std::vector<std::vector<int> > > node_edge_list(m_nPixels);

    // [edge_no, [first_node, second_node]]*
    std::vector<std::vector<int> > edge_list;

    extractEdgeList(&node_edge_list, &edge_list);

    // 0 indicates not used, 1 indicates edge is used by at least one tree.
    std::vector<int> global_used_edges(edge_list.size(), 0);

    // [number of trees][edges numbers]
    // ie [tree_no [edge_no]*]*
    std::vector<std::vector<int> > trees_edge_list;

    // sample usages :
    // extractPrimTreeEdges(0, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);
    // extractPrimTreeEdges(1, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);
    // extractDFSTreeMaxLimitNonRepeatableEdges(20, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);

    // you cannot use Prims right now because the later part of the code assumes that the first tree is a spanning tree
    // with no shared edges in the remaining trees.


    printf(" Creating Trees ... \n");
    while (areThereUnusedEdges(global_used_edges))
    {
      extractDFSTreeMaxLimitNonRepeatableEdges(edge_list.size()+1, node_edge_list, edge_list, &global_used_edges, &trees_edge_list);
    }



    // by using edge_list.size(), i am guaranteeing that the first tree is going to be a spanning tree.
    // this tree is used for getting the primal energies and the final solution.

    // after calling all trees dump edges for each tree for debugging purpose.
    // dumpTrees(trees_edge_list, edge_list);

    // this data structure is [graph_node_no][tree][0:tree_no, 1:tree_node_no]
    std::vector<std::vector<std::vector<int> > > node_in_trees(m_nPixels, (std::vector<std::vector<int> >)0);
    getWhatTreesNodeIsUsedIn(trees_edge_list, edge_list, &node_in_trees);

    for (unsigned int i=0; i < node_in_trees.size(); i++)
    {
      // make sure that each node is used at least in one tree.
      assert(node_in_trees[i].size() != 0);
    }


    const int num_trees = trees_edge_list.size();

    // It is important to note the node numbering used for trees.
    // ie what node numbers in the tree correspond to what node numbers
    // of the full graph
    // same for the edges.
    // hence we have mapping functions for both these.
    // the node number mapping is got from computeDataCostForTree
    std::vector<DataCost*> dcost_vec(num_trees, (DataCost*)0);
    SmoothnessCost* scost_vec = 0;
    std::vector<EnergyFunction*> energy_vec(num_trees, (EnergyFunction*)0);
    std::vector<MRF*> mrf_vec(num_trees, (MRF*)0);

    std::vector<MRF::CostVal*> D_a(trees_edge_list.size(), (MRF::CostVal*)0);
    for (unsigned int t=0; t<num_trees; t++)
    {
      // the number of nodes is edges+1
      D_a[t] = new MRF::CostVal[(trees_edge_list[t].size()+1)*m_nLabels];
    }

    std::vector<std::vector<int> > tree_node_to_graph_node(num_trees, (std::vector<int>)0);

    printf(" Computing DataCost for trees ... \n");
    computeDataCostForTrees(trees_edge_list, edge_list, node_in_trees, &D_a, &tree_node_to_graph_node, &dcost_vec);

    scost_vec = new SmoothnessCost((MRF::SmoothCostGeneralFn) m_smoothFn);

    // create mrfs of MaxProdBPTreeDD for each tree we have.
    for (int t=0; t < num_trees; t++)
    {
      energy_vec[t] = new EnergyFunction(dcost_vec[t], scost_vec);
      MaxProdBPTreeDD* mptd = new MaxProdBPTreeDD(trees_edge_list[t].size()+1, m_nLabels, energy_vec[t]);
      mptd->set_tree_node_to_graph_node(&tree_node_to_graph_node[t]);
      mrf_vec[t] = mptd;

      // set the edges for mrfs now
      for (int e=0; e < trees_edge_list[t].size(); e++)
      {
        int graph_edge_no = trees_edge_list[t][e];
        int graph_node1 = edge_list[graph_edge_no][0];
        int graph_node2 = edge_list[graph_edge_no][1];
        int tree_node1, tree_node2;
        getTreeNodeNo(node_in_trees[graph_node1], t, &tree_node1);
        getTreeNodeNo(node_in_trees[graph_node2], t, &tree_node2);
        mrf_vec[t]->setNeighbors(tree_node1, tree_node2, 1);
      }

      // Check validity of nodes connected in mrf.
      // will need to work on this.

      mrf_vec[t]->initialize();
      mrf_vec[t]->clearAnswer();

      // print output labels
      int num_tree_nodes = (trees_edge_list[t].size()+1);
      for (int i=0; i < num_tree_nodes; i++)
      {
        // printf("  %d", mrf_vec[t]->getLabel(i));
      }
      // printf("\n");

    }

    int changed = 1;
    int iter_no = 1;
    int max_iter = 300; // refers to infer from trees not gradient descent times
    // ie gradient descent times is 1 less than max_iter

    float gamma = 0.1;

    int gmismatchno = 0;
    std::vector<int> glabels(m_nPixels);


    MRF::EnergyVal e_trees_total_old = 0.0;
    MRF::EnergyVal e_trees_total = 0.0;
    MRF::EnergyVal best_primal_energy = 0.0;

    printf("\n*******  Iterating ****\n");


    while (changed == 1)
    {
      changed = 0;
      // when any label in any tree changes, set changed to 1.

      e_trees_total = 0.0;

      for (int t=0; t < num_trees; t++)
      {
        MRF::EnergyVal E;
        float time_t,tot_time_t;
        E = mrf_vec[t]->totalEnergy();
        tot_time_t= 0;
        mrf_vec[t]->optimize(1, time_t);
        E = mrf_vec[t]->totalEnergy();
        tot_time_t = tot_time_t + time_t ;

        e_trees_total += E;
      }

      // printf("**** Computing norm results from the trees ****\n");
      // I call the sum of result labels divide by number of trees that contain the node
      // as norm_results.
      std::vector<std::vector<float> > norm_results(m_nPixels, std::vector<float>(m_nLabels));
      // tree_count is stored in node_in_trees[i].size()

      // for each node
      for (unsigned int i=0; i<node_in_trees.size(); i++)
      {
        for (int t=0; t < node_in_trees[i].size(); t++)
        {
          int tree_no = node_in_trees[i][t][0];
          int tree_node_no = node_in_trees[i][t][1];
          int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);
          norm_results[i][label_no]++;
        }
      }

      // for each node
      for (unsigned int i=0; i<node_in_trees.size(); i++)
      {
        for (unsigned int d=0; d<m_nLabels; d++)
        {
          norm_results[i][d] /= node_in_trees[i].size();
        }
      }

      e_trees_total_old = e_trees_total;



      std::vector<int> mismatched_nodes; //(node_in_trees.size(), 0);
      int num_mismatch = 0;
      if (iter_no == 1)
      {
        changed = 1;
      } else
      {
        for (int ii=0; ii<m_nPixels; ii++)
        {
          setLabel(ii, -1); // Label allows -ve numbers.
        }
        // for each node

        for (unsigned int ii=0; ii<node_in_trees.size(); ii++)
        {
          int current_label = getLabel(ii);
          for (int t=0; t < node_in_trees[ii].size(); t++)
          {
            int tree_no = node_in_trees[ii][t][0];
            int tree_node_no = node_in_trees[ii][t][1];
            int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);
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
      // since iteration 1 starts with all tree nodes assigned to same value 0, mismatch will be 0
      // and so we should look at iteration 2
      // if (iter_no <= 2 || num_mismatch < gmismatchno)
      // {
      //  for (unsigned int iii=0; iii<m_nPixels; iii++)
      //  {
      //    glabels[iii] = getLabel(iii);
      //  }
      //  gmismatchno = num_mismatch;
      //}


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
        if (primal_energy < best_primal_energy || (fabs(primal_energy - best_primal_energy) < getToleranceParam() && gmismatchno > num_mismatch) )
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
        printf("** %d : mismatch = %d **  PE = %f  DE = %f ( %f secs ) ", iter_no, num_mismatch, primal_energy, e_trees_total, timess);
      // if (iter_no % 20 == 0)
        printf("\n");

      iter_no++;
      if (changed == 0 || iter_no > max_iter  || fabs(primal_energy - e_trees_total) < getToleranceParam())
        break;


      // computing the norm of the gradient
      float norm_grad = 0;

      for (unsigned int ii=0; ii<mismatched_nodes.size(); ii++)
       {
         int i = mismatched_nodes[ii];

         for (unsigned int t=0; t < node_in_trees[i].size(); t++)
         {
           int tree_no = node_in_trees[i][t][0];
           int tree_node_no = node_in_trees[i][t][1];
           int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);

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

      std::vector<int> trees_needed_to_be_updated(num_trees, 0);

      //for (unsigned int i=0; i<node_in_trees.size(); i++)
      // clock_t startt = clock();
      for (unsigned int ii=0; ii<mismatched_nodes.size(); ii++)
      {
        int i = mismatched_nodes[ii];
        // printf("**** i = %d \n", i);

        for (unsigned int t=0; t < node_in_trees[i].size(); t++)
        {
          int tree_no = node_in_trees[i][t][0];
          trees_needed_to_be_updated[tree_no] = 1;
          int tree_node_no = node_in_trees[i][t][1];
          int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);
          // This x_i_treeno(label_no) = 1 all else x_i_treeno are 0
          // norm_results[i][label_no]++;
          // MaxProdBPTreeDD* mrf_temp = (MaxProdBPTreeDD*)mrf_vec[tree_no];
          int idx;

          for (unsigned int d=0; d<m_nLabels; d++)
          {
            idx = tree_node_no*m_nLabels + d;
            float theta_p_d = mrf_vec[tree_no]->getEnergyFunction()->m_dataCost->getDataCostRaw(idx);
            float gradient = -norm_results[i][d];
            if (d == label_no)
              gradient += 1;
            // else x_p_d is 0

            gradient *= gamma*(best_primal_energy - e_trees_total)/norm_grad;
            mrf_vec[tree_no]->getEnergyFunction()->m_dataCost->updateDataCostRaw(idx, theta_p_d + gradient);
            // mrf_temp->updateData(mrf_vec[tree_no]->getEnergyFunction()->m_dataCost->getDataCostArray());
          }

        }
      }

      for (int i=0; i < num_trees; i++)
      {
        if (trees_needed_to_be_updated[i] == 1)
        {
          MaxProdBPTreeDD* mrf_temp = (MaxProdBPTreeDD*)mrf_vec[i];
          mrf_temp->updateData(mrf_vec[i]->getEnergyFunction()->m_dataCost->getDataCostArray());
        }
      }

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


    for (unsigned int t=0; t<num_trees; t++)
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


// Use this when you are sure first tree is a spanning tree
/*
for (unsigned int i=0; i<node_in_trees.size(); i++)
{
  int current_label = getLabel(i);
  int flag_tree_assigned = 0;
  for (int t=0; t < node_in_trees[i].size(); t++)
  {
    int tree_no =  node_in_trees[i][t][0];
    if (tree_no != 0)
      continue;
    flag_tree_assigned = 1;
    int tree_node_no = node_in_trees[i][t][1];
    int label_no = mrf_vec[tree_no]->getLabel(tree_node_no);
    // since we are only looking at tree 0,
    current_label = label_no;
  }
  if (flag_tree_assigned == 0)
  {
    printf(" Tree 0 not spanning error \n");
    exit(1);
  }
  if (current_label == -1)
    printf(" Serious representation error \n");
  else
    setLabel(i, current_label);
}
*/



