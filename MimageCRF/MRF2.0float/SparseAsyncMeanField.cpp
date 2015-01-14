#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ArrayMath.h"
#include "SparseAsyncMeanField.h"

#include <time.h>

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeMessages() -> 
 *  computeUpDownMessages(),computeLeftRightMessages() -> 
 *  computeMessage()
 */

#define PIXEL_IND(r,c) (r)*m_width+(c)

SparseAsyncMeanField::SparseAsyncMeanField(int width, int height, int nLabels, 
                             EnergyFunction *eng) :
  SparseSyncMeanField(width,height,nLabels,eng)
{ }

SparseAsyncMeanField::SparseAsyncMeanField(int nPixels, int nLabels,
                                         EnergyFunction *eng) :
  SparseSyncMeanField(nPixels,nLabels,eng)
{  
  *m_log << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

SparseAsyncMeanField::~SparseAsyncMeanField() {}

void SparseAsyncMeanField::initializeAlg(InitBelType init)
{

  switch(init)
    {
    case UNIFORM:
      SparseSyncMeanField::initializeAlg ();
      break;
    case LOCAL:
      SyncMeanField::initializeAlg(init);
      //updateSparsity();
      break;
    case LABELPOINT:
      {
        FloatType epsilon = 1e-8;
        FloatType logProb0 = log(epsilon) - log((FloatType)m_nLabels-1);
        FloatType logProb1 = log(1-epsilon);
        
        FloatType prob0 = epsilon/m_nLabels;
        FloatType prob1 = 1-epsilon;

        for (int i=0 ; i<m_nPixels ; i++)
          {
            for (int b=0 ; b<m_nLabels ; b++)
              {
                m_nodeBeliefs[i][b] = logProb0;
                m_expNodeBeliefs[i][b] = prob0;
                m_spExpNodeBeliefs[i][b] = 0;
                m_nodeSparsity[i][b] = false;
              }

            m_nodeBeliefs[i][m_labels[i]] = logProb1;
            m_expNodeBeliefs[i][m_labels[i]] = prob1;

            m_spExpNodeBeliefs[i][m_labels[i]] = 1.0;
            m_nodeSparsity[i][m_labels[i]] = true;
          }
        break;
      }
    }

  initAnswer();  

}

/** Run the algorithm (asynchronous sparse mean field)
 */     

void SparseAsyncMeanField::optimizeAlg (int nIterations)
{

  FloatType diff, energy = 0;

  FloatType energy0;

  if (m_beliefDiffFun==ENERGY)
    energy0 = sparseFreeEnergy();
  else
    energy0 = 0;

  FloatType energyOld = energy0;
  FloatType time = 0;
  
  //*m_log << "Initial free energy " << energy0 << "\n";

  //// Print initial free energy and info ////
  *m_log << -1 << "\tNaN\t";

  m_log->precision(1);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << energy0 <<"\t";

  m_log->precision(2);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << time << endl;

  if (m_imWriter!=0)
    {
      updateAnswer();
      m_imWriter->write(m_labels);
    }

  // Print marginal entropies
//    *m_log <<"H ";
//    for (int p=0 ; p<m_nPixels ; p++)
//      *m_log << getNodeEntropy(p) <<" ";
//    *m_log << endl;

  // Print sparsity
//    *m_log <<"S ";
//    for (int p=0 ; p<m_nPixels ; p++)
//      *m_log << getNumNonSparseStates(p) <<" ";
//    *m_log << endl;


  clock_t start, end;

  for (int i=0 ; i<nIterations ; i++)
    {

      start = clock();

      diff = computeUpdateBeliefs();
  
      // Print marginal entropies
//        *m_log <<"H ";
//        for (int p=0 ; p<m_nPixels ; p++)
//          *m_log << getNodeEntropy(p) <<" ";
//        *m_log << endl;


       // Print sparsity
//        *m_log <<"S ";
//        for (int p=0 ; p<m_nPixels ; p++)
//          *m_log << getNumNonSparseStates(p) <<" ";
//        *m_log << endl;
       
      if (m_beliefDiffFun==ENERGY)
        {
          energy = sparseFreeEnergy();
          diff = 100*(energyOld-energy)/energyOld;
          energyOld = energy;
        }
      else
        std::cout<<"it ain't energy!\n";

      end = clock();

      time += (float) (((double)(end-start)) / CLOCKS_PER_SEC);
      
//       *m_log << "SparseAsyncMeanField: Iter = " << i << "\tDiff = " << diff <<
//         "\tFree Energy = " << energy << "\n";

      //// Print iteration free energy and info ////
      *m_log << i << "\t";

      m_log->precision(5);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << diff << "\t";

      m_log->precision(1);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << energy << "\t";

      m_log->precision(2);
      m_log->setf(std::ios::fixed,std::ios::floatfield);
      *m_log << time << endl;
      ////

//       *m_log<< "Average sparsity = " << 
//         ((float)getNumNonSparseStates())/((float)m_nPixels*m_nLabels)<<'\n';


//       printHistSparsity(10);
      
      if (m_imWriter!=0)
        {
          updateAnswer();
          m_imWriter->write(m_labels);
        }

      if (diff < m_beliefDiffTol)
        {
//           *m_log << "SparseAsyncMeanField converged in " << i << 
//             " iterations.\nFree energy reduction " << 
//             100*(energy0-energy)/energy0 << "%.\n";

          break;
        }
    }

  updateAnswer();

#ifndef NDEBUG
  //  printAll();
  //  printAllBeliefs();
#endif


}


/** Asynchronously compute all new beliefs from the current beliefs.
 */
FloatType SparseAsyncMeanField::computeUpdateBeliefs()
{
  
  FloatType diff,diff0;

  diff = 0;
  int pixel;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        computeBelief (r,c);

        pixel = PIXEL_IND(r,c);

        diff0 = beliefDifference (pixel);

        diff = combineDifference (diff,diff0);

        updateBelief (pixel);
        updateSparsity (pixel);

      }

  return diff;
}


FloatType SparseAsyncMeanField::sparseFreeEnergy()
{

  register FloatType energy = 0;

  int pixel,pixel1;
  
  FloatType *belief, *belief1;
  bool *sp, *sp1;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        pixel = PIXEL_IND(r,c);

        // Data terms
        belief = m_spExpNodeBeliefs[pixel];
        sp = m_nodeSparsity[pixel];

        for (int b=0 ; b<m_nLabels ; b++ )
          if (sp[b])
            energy += belief[b] * getDataCost(pixel,b);

        // Smoothness terms

        if (r<m_height-1)
          { // Down
            pixel1 = PIXEL_IND(r+1,c);

            belief = m_spExpNodeBeliefs[pixel];
            belief1 = m_spExpNodeBeliefs[pixel1];
            
            sp = m_nodeSparsity[pixel];
            sp1 = m_nodeSparsity[pixel1];

            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                if (sp[b] && sp1[b1])
                  energy += belief[b] * belief1[b1] * 
                    getSmoothnessCost(pixel,pixel1,b,b1);
            
          }

        if (c<m_width-1)
          { // Right
            pixel1 = PIXEL_IND(r,c+1);

            belief = m_spExpNodeBeliefs[pixel];
            belief1 = m_spExpNodeBeliefs[pixel1];
            
            sp = m_nodeSparsity[pixel];
            sp1 = m_nodeSparsity[pixel1];
            
            for (int b=0 ; b<m_nLabels ; b++ )
              for (int b1=0 ; b1<m_nLabels ; b1++)
                if (sp[b] && sp1[b1])
                  energy += belief[b] * belief1[b1]  * 
                    getSmoothnessCost(pixel,pixel1,b,b1);
            
          }
            

      }

  return energy - getSparseTotalNodeEntropy();
    

  
  
}

FloatType SparseAsyncMeanField::getSparseTotalNodeEntropy()
{

  register FloatType entropy = 0;

  FloatType** beliefs = m_spExpNodeBeliefs;
  FloatType* bel;

  bool** sparsity = m_nodeSparsity;
  bool* sp;
  for (int i=0 ; i<m_nPixels ; i++, beliefs++, sparsity++)
    {
      bel = *beliefs;
      sp = *sparsity;

      for (int b=0 ; b<m_nLabels ; b++, bel++, sp++)
        if (*sp)
          entropy +=  (*bel) * log(*bel);
    }
  
  return -entropy;
}
