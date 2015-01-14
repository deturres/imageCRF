#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ArrayMath.h"
#include "MeanField.h"

using std::endl;

/** 
 * You should be able to trace the general flow through:
 *  optimizeAlg() -> 
 *  computeUpdateBeliefs() -> 
 *  computeBelief()
 */



MeanField::MeanField(int width, int height, int nLabels, 
                             EnergyFunction *eng) :
  MRFEnergy(width,height,nLabels,eng)
{
  initializeMeanField();
}

MeanField::MeanField(int nPixels, int nLabels,EnergyFunction *eng) :
  MRFEnergy(nPixels,nLabels,eng)
{  
  std::cerr << "Non grid graphs are not supported";
  assert(m_grid_graph);
}

MeanField::~MeanField() 
{
  for (int i=0 ; i<m_nPixels ; i++)
    {
      delete[] m_nodeBeliefs[i];
      delete[] m_expNodeBeliefs[i];
      delete[] m_newNodeBeliefs[i];
    }
  
  delete[] m_nodeBeliefs;
  delete[] m_expNodeBeliefs;
  delete[] m_newNodeBeliefs;
}


/** Common initialization routine */
void MeanField::initializeMeanField()
{

  m_beliefDamper = 1.0;
  m_beliefDiffFun = KLD;
  m_beliefDiffTol = .001;

  // Belief storage
  m_nodeBeliefs = new FloatType*[m_nPixels];
  m_newNodeBeliefs = new FloatType*[m_nPixels];
  m_expNodeBeliefs = new FloatType*[m_nPixels];

  //
  // Create beliefs
  //
  
  for (int i=0 ; i<m_nPixels ; i++)
    {
      m_nodeBeliefs[i] = new FloatType[m_nLabels];
      m_newNodeBeliefs[i] = new FloatType[m_nLabels];
      m_expNodeBeliefs[i] = new FloatType[m_nLabels];
    }

}

/** Initialization from parent classes.
 *
 * Start beliefs at uniform.
 */
void MeanField::initializeAlg()
{
  initializeAlg(UNIFORM);
}

/** Initialize beliefs */
void MeanField::initializeAlg(InitBelType init)
{

  switch (init)
    {
    case UNIFORM:
      {
        FloatType logZ = -log(m_nLabels);
        FloatType Z = 1.0/((FloatType)m_nLabels);
        
        for (int i=0 ; i<m_nPixels ; i++)
          for (int b=0 ; b<m_nLabels ; b++)
            {
              m_nodeBeliefs[i][b] = logZ;
              m_expNodeBeliefs[i][b] = Z;
            }
        break;
      }
    case LOCAL:
      {
        for (int i=0 ; i<m_nPixels ; i++)
          {
            copyDataCost (i,m_nodeBeliefs[i]);
            logNormalize(m_nodeBeliefs[i]);

            for (int b=0 ; b<m_nLabels ; b++)
              m_expNodeBeliefs[i][b] = exp(m_nodeBeliefs[i][b]);
      
          }
        break;
      }
    case LABELPOINT:
      {
        //FloatType negInf = - std::numeric_limits<FloatType>::infinity();
        FloatType epsilon = 1e-8;
        FloatType logProb0 = log(epsilon) - log(m_nLabels-1);
        FloatType logProb1 = log(1-epsilon);

        FloatType prob0 = epsilon/m_nLabels;
        FloatType prob1 = 1-epsilon;

        for (int i=0 ; i<m_nPixels ; i++)
          {
            for (int b=0 ; b<m_nLabels ; b++)
              {
                m_nodeBeliefs[i][b] = logProb0;
                m_expNodeBeliefs[i][b] = prob0;
              }

            m_nodeBeliefs[i][m_labels[i]] = logProb1;
            m_expNodeBeliefs[i][m_labels[i]] = prob1;
          }
        break;
      }
    }
  
  initAnswer();
}

/** Run the algorithm ( plain-vanilla mean field)
 */     

void MeanField::optimizeAlg (int nIterations)
{

  FloatType diff,energy = 0;
  
  FloatType energy0;

  if (m_beliefDiffFun==ENERGY)
    energy0 = freeEnergy();
  else
    energy0 = 0;
  
  FloatType energyOld = energy0;

  FloatType time = 0;
  
  // Log format is (tab separated): iteration energy time

  *m_log << -1 << "\tNaN\t";

  m_log->precision(1);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << energy0 <<"\t";

  m_log->precision(2);
  m_log->setf(std::ios::fixed,std::ios::floatfield);
  *m_log << time << endl;

  // Write intermediate image
  if (m_imWriter!=0)
    {
      updateAnswer();
      m_imWriter->write(m_labels);
    }

  
  clock_t start, end;
  
  for (int i=0 ; i<nIterations ; i++)
    {

      start = clock();

      diff = computeUpdateBeliefs();

      if (m_beliefDiffFun==ENERGY)
        {
          energy = freeEnergy();
          diff = 100*(energyOld-energy)/energyOld;
          energyOld = energy;
        }

      end  = clock();

      time += (float) (((double)(end-start)) / CLOCKS_PER_SEC);
      
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

      if (m_imWriter!=0)
        {
          updateAnswer();
          m_imWriter->write(m_labels);
        }
      
      if (diff < m_beliefDiffTol)
          break;
    }

  updateAnswer();

}

/** Asynchronously compute all new beliefs from the current beliefs.
 */
FloatType MeanField::computeUpdateBeliefs()
{
  
  FloatType diff,diff0;

  diff = 0;
  int pixel;

  for (int r=0 ; r<m_height ; r++)
    for (int c=0 ; c<m_width ; c++)
      {
        computeBelief (r,c);

        pixel = PIXEL_IND(r,c);

        diff0 = beliefDifference(pixel);

        diff = combineDifference(diff,diff0);

        updateBelief(pixel);
        

      }

  return diff;
}

/** Compute a new belief at a pixel from the current beliefs.
 */
void MeanField::computeBelief(int r, int c) 
{
  int pixC = pixelIndex(r,c);
  int pixL = pixelIndex(r,c-1);
  int pixR = pixelIndex(r,c+1);
  int pixU = pixelIndex(r-1,c);
  int pixD = pixelIndex(r+1,c);
        
  // Get the belief buffer
  FloatType* belief = m_newNodeBeliefs[pixC];
  
  // Start with -D(p)
  copyDataCost(pixC, belief);
  
  if (c>0)
    // Left neighbor
    addMeanSmoothnessCost(pixC, pixL, belief);
  
  if (c<m_width-1)
    // Right neighbor
    addMeanSmoothnessCost(pixC, pixR, belief);
  
  if (r>0)
    // Up neighbor
    addMeanSmoothnessCost(pixC, pixU, belief);
  
  if (r<m_height-1)
    // Down neighbor
    addMeanSmoothnessCost(pixC, pixD, belief);
  
  logNormalize(belief);
  
}


/** Add the mean smoothness cost to a vector using the current beliefs
 *
 * Parameters:
 *  - dst is a buffer over which fixedPixel indexes
 *
 * Output: M[i] = - sum_j b[j] * V(i,j)
 *
 */
void SyncMeanField::addMeanSmoothnessCost(const int fixedPixel, 
                                          const int varPixel, FloatType* dst)
{
  
  FloatType* belief0 = m_expNodeBeliefs[varPixel];
  FloatType* belief;

  register FloatType mean;

  for (int i=0 ; i<m_nLabels ; i++, dst++)
    {
      belief = belief0;
      
      mean = 0;

      for (int j=0 ; j<m_nLabels ; j++, belief++)
        {
          mean += (*belief) * getSmoothnessCost(fixedPixel, varPixel, i, j);
        }

      *dst -= mean;
    }
}
