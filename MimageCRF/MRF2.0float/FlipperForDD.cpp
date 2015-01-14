/*
 * FlipperForDD.cpp
 *
 *  Created on: May 18, 2012
 *      Author: bhole
 */

#include <time.h>
#include <limits.h>
#include <math.h>
#include <queue>

#include "FlipperForDD.h"
#include "GCoptimization.h"



FlipperForDD::FlipperForDD(int width, int height, int nLabels,EnergyFunction *eng):MRFWrapper(width,height,nLabels,eng)
{
}

FlipperForDD::FlipperForDD(int nPixels, int nLabels,EnergyFunction *eng):MRFWrapper(nPixels,nLabels,eng)
{
}

FlipperForDD::~FlipperForDD()
{
  if (m_grid_graph)
  {
    printf(" Flipper does not handle grid graphs interface \n");
    exit(1);
  } else {
    // for non-grid graphs


    // for non-grid graphs

    // exact inference
    // create flipper graph and run graphcuts

    if (m_nLabels != 2)
    {
      printf(" Current algorithm does not work for other than binary case. \n");
      return;
    }

    delete energy_flipper;
    delete mrf_flipper;
    delete sscost;
    delete dcost_flipper;
    delete[] D_a;

  }


}


void FlipperForDD::setFlipperFunction(MRF::SmoothCostGeneralFn fnCost_flipper)
{
  m_flippersmoothFn = fnCost_flipper;
}


void FlipperForDD::setflippedinfo(std::vector<int> **flippedinfo_f)
{
  flippedinfo_ptop = flippedinfo_f;
}


void FlipperForDD::extractEdgeList(std::vector<std::vector<std::vector<int> > > *node_edge_list,
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

int FlipperForDD::checkValidFlipperGraph()
{
  return 0;
}


// mark them in node_in_subgraphs
// now we always need to check the edges to know what nodes are flipped.
void FlipperForDD::markFlippedNodes(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<int> &edge_list_type,
    const std::vector<std::vector<int> > &edge_list,
    std::vector<int> *nodes_flipped)
{

  std::queue<int> flippedlist;
  std::vector<int> checked(m_nPixels, 0);
  // because there can be forests of graphs, we need to go thru all nodes

  // Let the first pixel of each graph in a forest always be non-flipped
  // so i don't need a seeding strategy.

  for (int i=0; i<m_nPixels; i++)
  {
    if (checked[i] == 0 && node_edge_list[i].size() > 0) // this allows start of a new forest
    {
      // i is the first pixel of a new graph in a forest, so it is non-flipped
      checked[i] = 1;
      for (int nn=0; nn < node_edge_list[i].size(); nn++)
      {
        int neighbor = node_edge_list[i][nn][0];
        int edge = node_edge_list[i][nn][1];
        if (checked[neighbor] == 0)
        {
          flippedlist.push(neighbor);
          checked[neighbor] = 1;

          // for the first node in the graph in the forest, we can only check
          // if the edge is nonsubmodular since we know the first node is not flipped
          // but for completeness we check both conditions
          if (edge_list_type[edge] == 1 && (*nodes_flipped)[i] == 0) // nonsubmodular and parent not flipped
          {
            (*nodes_flipped)[neighbor] = 1;
          } else if (edge_list_type[edge] == 0 && (*nodes_flipped)[i] == 1) // submodular but parent is flipped
          {
            (*nodes_flipped)[neighbor] = 1;
          }

        } else {
          printf(" This was a new tree, so this is not possible \n");
          exit(1);
        }

      }

      while (!flippedlist.empty())
      {
        // nodep stands for parent node
        int nodep = flippedlist.front();
        flippedlist.pop();

        for (int nn=0; nn < node_edge_list[nodep].size(); nn++)
        { // go through all neighbors
          int neighbor = node_edge_list[nodep][nn][0];
          int edge = node_edge_list[nodep][nn][1];

          if (checked[neighbor] == 1)
          {
            // check consistency of flipper graph is maintained otherwise,
            // throw an error
            int type_nodep = (*nodes_flipped)[nodep];
            int type_nodec = (*nodes_flipped)[neighbor];
            if (edge_list_type[edge] == 1 && type_nodep == 0 && type_nodec == 1) {}
            else if (edge_list_type[edge] == 1 && type_nodep == 1 && type_nodec == 0) {}
            else if (edge_list_type[edge] == 0 && type_nodep == 0 && type_nodec == 0) {}
            else if (edge_list_type[edge] == 0 && type_nodep == 1 && type_nodec == 1) {}
            else {
              printf(" This is not a flipper graph \n");
              exit(1);
            }

          } else if (checked[neighbor] == 0)
          {
            flippedlist.push(neighbor);
            checked[neighbor] = 1;
            int type_nodep = (*nodes_flipped)[nodep];

            if (edge_list_type[edge] == 1 && type_nodep == 0) // nonsubmodular and parent not flipped
            {
              (*nodes_flipped)[neighbor] = 1;
            } else if (edge_list_type[edge] == 0 && type_nodep == 1) // submodular but parent is flipped
            {
              (*nodes_flipped)[neighbor] = 1;
            }
          }
        }
      }
    }
  }
}


// the node number order is maintained.
void FlipperForDD::computeDataCostForFlipper(
    const std::vector<std::vector<std::vector<int> > > &node_edge_list,
    const std::vector<int> &edge_list_type,
    const std::vector<std::vector<int> > &edge_list,
    const std::vector<int> &nodes_flipped,
    MRF::CostVal* D_a)
{

  assert(m_nLabels == 2);

  int dsiIndex = 0;

  for (unsigned int i=0; i < m_nPixels; i++)
  {
    if (nodes_flipped[i] == 0) // not flipped
    {
      for (int d = 0; d < m_nLabels; d++)
      {
        MRF::CostVal* ptr = D_a + dsiIndex;
        *ptr = nodeArray[i].localEv[d];
        dsiIndex++;
      }
    } else if (nodes_flipped[i] == 1) // flipped
    {
      for (int d = m_nLabels-1; d >= 0; d--)
      {
        MRF::CostVal* ptr = D_a + dsiIndex;
        *ptr = nodeArray[i].localEv[d];
        dsiIndex++;
      }

    }
  }

}



void FlipperForDD::updateFlipperDataStructures(const int idx, const int pixelno, const int d,
    const float value)
{

  // no need to consider flipping for Flipper, since this is the outer Flipper data structure
  getEnergyFunction()->m_dataCost->updateDataCostRaw(idx, value);

  // i probably also need to update : nodeArray[i].localEv[d]
  // note that this doesn't seem to have an effect but good to do so all structures
  // are kept updated
  // flipping does not seem to be correct here,
  // because in function above computeDataCostForFlipper, we copy from [d]th location
  // in nodeArray to D_a flipped ie d if not flipped, 1-d if flipped
  nodeArray[pixelno].localEv[d] = value;

  // need to update the datastructure used by Expansion() because setUp Flipper is called
  // only once, and since we update datacost terms in between, we need to go directly inside
  // note that Expansion will have flipped nodes and so we need to take care about flipping
  assert(m_nLabels == 2);
  int idx2 = idx;
  if (nodes_flipped[pixelno] == 0) {} // not flipped, do nothing
  else if (nodes_flipped[pixelno] == 1)
  {
    if (d == 0)
      idx2 ++; // make it di = 1;
    else if (d==1)
      idx2--;  // make it di = 0;
  }
  mrf_flipper->getEnergyFunction()->m_dataCost->updateDataCostRaw(idx2, value);

}


void FlipperForDD::setUpFlipperForDD()
{

  if (m_grid_graph)
  {
    printf(" Flipper does not handle grid graphs interface \n");
    exit(1);
  } else {
    // for non-grid graphs


    // for non-grid graphs

    // exact inference
    // create flipper graph and run graphcuts

    if (m_nLabels != 2)
    {
      printf(" Current algorithm does not work for other than binary case. \n");
      return;
    }


    node_edge_list.resize(m_nPixels);

    extractEdgeList(&node_edge_list, &edge_list);

    edge_list_type.resize(edge_list.size(),0);

    // this notes what edges are submodular and what are nonsubmodular
    extractEdgeType(edge_list, &edge_list_type);

    // will at some point need to check that the graph is a valid Flipper graph
    // checkValidFlipperGraph();

    nodes_flipped.resize(m_nPixels, 0);
    markFlippedNodes(node_edge_list, edge_list_type, edge_list, &nodes_flipped);

    dcost_flipper = NULL;
    energy_flipper = NULL;
    mrf_flipper = NULL;
    D_a = new MRF::CostVal[m_nPixels*m_nLabels];

    //printf(" Computing DataCost for Flipper graph ... \n");
    computeDataCostForFlipper(node_edge_list, edge_list_type, edge_list, nodes_flipped, D_a);
    dcost_flipper = new DataCost (D_a);

    sscost = new SmoothnessCost((MRF::SmoothCostGeneralFn) m_flippersmoothFn);
    energy_flipper = new EnergyFunction(dcost_flipper, sscost);
    mrf_flipper = new Expansion(m_nPixels, m_nLabels, energy_flipper);

    // set the edges for mrfs now
    for (int e=0; e < edge_list.size(); e++)
    {
      int graph_node1 = edge_list[e][0];
      int graph_node2 = edge_list[e][1];
      mrf_flipper->setNeighbors(graph_node1, graph_node2, 1);
    }

    mrf_flipper->initialize();
    mrf_flipper->clearAnswer();

  }



}



void FlipperForDD::optimizeAlg(int nIterations)
{

  if (m_grid_graph)
  {
    printf(" Flipper does not handle grid graphs interface \n");
    exit(1);
  } else {
    // for non-grid graphs

    // exact inference
    // create flipper graph and run graphcuts

    if (m_nLabels != 2)
    {
      printf(" Current algorithm does not work for other than binary case. \n");
      return;
    }


    // stop timer
    clock_t finishss = clock();
    float timess = (float) (((double)(finishss - start)) / CLOCKS_PER_SEC);
    // printf(" Time for creating and initializing flipper graphs = %f secs\n", timess);

    // VERY IMPORTANT
    // Every time you call mrf_vec[t]->optimize(1, time_t); or mrf_vec[t]->totalEnergy();
    // you need to make sure that
    // the correct flipstructure and subgraph_node_to_graph_node structure are activated
    // for the external function

    MRF::EnergyVal E;
    float time_t=0;
    *flippedinfo_ptop = &nodes_flipped;

    //E = mrf_flipper->totalEnergy();
    mrf_flipper->optimize(1, time_t);
    E = mrf_flipper->totalEnergy();

    for (int i=0; i<m_nPixels; i++)
    {
      int label_no = mrf_flipper->getLabel(i);
      if (nodes_flipped[i] == 1)
        label_no = 1 - label_no;
      setLabel(i, label_no);
    }

    finishss = clock();
    timess = (float) (((double)(finishss - start)) / CLOCKS_PER_SEC);

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

  }

}
