#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "ArrayMath.h"
#include "MeanField.h"

using std::endl;

/**
 * You should be able to trace the general flow through:
 *  optimizeAlg() ->
 *  computeUpdateBeliefs() ->
 *  computeBelief(), updateBelief()
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
	initializeMeanField();
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
    case CURRENT:
      break;
    case RESET:
      // Don't think there's anything to do here
      break;
    case UNIFORM:
      {
        FloatType logZ = -log((FloatType)m_nLabels);
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
        FloatType logProb0 = log(epsilon) - log((FloatType)m_nLabels-1);
        FloatType logProb1 = log(1-epsilon);

        FloatType prob0 = epsilon/(m_nLabels-1);
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


  if (m_grid_graph) {

  } else { //general graph

      m_maxnNB=0;

  	// traverse Neighbors list
  	for (int m=0; m < m_nPixels; m++) {

  		int numNB=traverseNeighbors((int)m);

  		m_nNeighbors[m] = numNB;

  		if (m_maxnNB<numNB) {
  			m_maxnNB=numNB;
  		}
  	}
  }

}

/** Run the algorithm ( plain-vanilla mean field)
 */

void MeanField::optimizeAlg (int nIterations)
{
  FloatType diff = 0,energy = 0;

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
          //m_imWriter->write(m_labels);
        }

      if (isnan(diff)) {
      } else {
		  if (diff < m_beliefDiffTol)
			  break;
      }
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

  if (m_grid_graph) {

	  for (int r=0 ; r<m_height ; r++)
		for (int c=0 ; c<m_width ; c++)
		  {
			computeBelief (r,c);

			pixel = pixelIndex(r,c);

			diff0 = beliefDifference(pixel);

			diff = combineDifference(diff,diff0);

			updateBelief(pixel);
		  }

  } else { //general graphs

	  for (pixel=0; pixel<m_nPixels; pixel++)
	  {
		  computeBelief (pixel);

		  diff0 = beliefDifference(pixel);

		  diff = combineDifference(diff,diff0);

		  updateBelief(pixel);
	  }

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

void MeanField::computeBelief(int pixel)
{
  // Get the belief buffer
  FloatType* belief = m_newNodeBeliefs[pixel];

  // Start with -D(p)
  copyDataCost(pixel, belief);

  for (int p=0; p<m_nNeighbors[pixel]; p++)
  {
	int pixel1=findNeighborWithIdx(pixel, p);
	if (pixel1==-1) {
		printf("Error\n");
		exit(1);
	}

	addMeanSmoothnessCost(pixel, pixel1, belief);
  }

  logNormalize(belief);
}








/** Update a local variational distribution based on the neighbors */
void MeanField::updateBelief(int pixel)
{
  //
  // Update log-space beliefs
  //
  if (m_beliefDamper==1)
    {
      // Update by swapping pointers
      FloatType* tmp = m_nodeBeliefs[pixel];

      m_nodeBeliefs[pixel] = m_newNodeBeliefs[pixel];
      m_newNodeBeliefs[pixel] = tmp;
    }
  else
    {
      // Update by convex combination
      FloatType *belief = m_nodeBeliefs[pixel];
      FloatType *nbelief = m_newNodeBeliefs[pixel];

      convexLogAdd (m_beliefDamper, belief, nbelief);
    }

  //
  // Update/calculate exponentiated beliefs
  //
  FloatType *ebelief = m_expNodeBeliefs[pixel];
  FloatType *belief = m_nodeBeliefs[pixel];

  for (int b=0 ; b<m_nLabels ; b++, ebelief++, belief++)
    *ebelief = exp(*belief);

}

/** Add the mean smoothness cost to a vector using the current beliefs
 *
 * Parameters:
 *  - dst is a buffer over which fixedPixel indexes
 *
 * Output: M[i] = - sum_j b[j] * V(i,j)
 *
 */
void MeanField::addMeanSmoothnessCost(const int fixedPixel,
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


/* Calculate the difference between the current and updated beliefs.
 *
 */

FloatType MeanField::beliefDifference (int pixel)
{
  FloatType* pb = m_nodeBeliefs[pixel];
  FloatType* pe = m_expNodeBeliefs[pixel];
  FloatType* pn = m_newNodeBeliefs[pixel];


  register FloatType diff0, diff = 0;

  switch (m_beliefDiffFun)
    {
    case NOP:
      break;
    case L1: // |p-q|

      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        diff += fabs( *pb - *pn );

      break;
    case L2: // (p-q)^2

      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        {
          diff0 = fabs( *pb - *pn );
          diff += diff0*diff0;
        }

        break;

    case MAX: // max_i |p_i - q_i|

      diff = fabs ( *pb - *pn );

      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++)
        {
          diff0 = fabs( *pb - *pn );
          if (diff0>diff)
            diff = diff0;
        }
      break;

    case KLD: // sum_i p_i * log(p_i/q_i)

      for (int b=0 ; b<m_nLabels ; b++, pb++, pn++, pe++)
        diff += (*pe) * ( *pb - *pn );

      break;

    case ENERGY:
      // Cannot calculated incrementally
      diff = 0;
      break;

    default:
      std::cerr<<"MeanField::beliefDifference " << m_beliefDiffFun << " not implemented yet. Ignoring.\n";
    }
  return diff;
}

/* Helper function for merging the message differences from two directions
 */
FloatType MeanField::combineDifference(FloatType diff1, FloatType diff2)
{
  FloatType diff = 0;

  switch (m_beliefDiffFun)
    {
    case NOP:
    case ENERGY:
      break;
    case L1:
    case L2:
    case KLD:
      diff = diff1+diff2;
      break;
    case MAX:
      diff = (diff2>diff1 ? diff2 : diff1);
      break;

    }

  return diff;
}

/** Convex log update, then difference
 * a - log((1-alpha)*exp(a)+alpha*exp(b))
 */
inline FloatType MeanField::convexLogDiff(const FloatType a, const FloatType b)
{
  if (m_beliefDamper==1)
    return a-b;
  else
    return a - log((1-m_beliefDamper)*exp(a)+m_beliefDamper*exp(b));
}

 /** Sets answer  to current Maximum Posterior Marginal (MPM) estimate
 */
void MeanField::updateAnswer()
{
   Label* label = m_labels;

   FloatType** belief = m_nodeBeliefs;

   for (int i=0; i<m_nPixels ; i++, label++, belief++)
   {
     *label = argMax (*belief);
     // (*label > 0 && i>= 30000) {
     /*
if(i>=30000) {
     std::cout << " Pix No: " << i << "  " ;
	 for (int mm=0; mm<m_nLabels ; mm++)
	 {
		 std::cout << (*belief)[mm] << "  ";
	 }
	 std::cout << std::endl;
     }  */
   }
}

/** Sets answer using only data cost */
void MeanField::initAnswer()
{
  Label* label = m_labels;

  FloatType *cost = new FloatType[m_nLabels];
  assert(cost!=0);
  for (int i=0 ; i<m_nPixels ; i++, label++)
    {
      copyDataCost(i, cost);
      *label = argMax (cost);
    }
  delete[] cost;
}

/** Calculate the current approximate marginal for a pair of neighboring pixels
 *  and write it to the destination array provided (must be of length
 * m_nlabels^2).
 */
void MeanField::getEdgeBelief(const int pixel1, const int pixel2,
                                FloatType* dst)
{
  FloatType* bel1 = m_nodeBeliefs[pixel1];
  FloatType* bel2 = m_nodeBeliefs[pixel2];

  //
  // Pixel 2
  //

  // Repeat the one dimensional belief over the ROWS of the destination

  FloatType *pd = dst;

  for (int b1=0 ; b1<m_nLabels ; b1++, pd+=m_nLabels)
    copyTo(pd, bel2);

  //
  // Pixel 1
  //

  // Repeat (add) the one dimensional belief over the COLS of the destination

  pd = dst;

  for (int b1=0 ; b1<m_nLabels ; b1++, bel1++)
    for (int b2=0 ; b2<m_nLabels ; b2++, pd++)
      *pd += *bel1;
}

/** Calculate the log probability of a labeling from the current beliefs
 */
FloatType MeanField::getLogProb(Label* labels)
{

  register FloatType logprob = 0;

  FloatType** belief = m_nodeBeliefs;

  for (int i=0 ; i<m_nPixels ; i++, belief++, labels++)
    {
      logprob += *((*belief) + (*labels));
    }

  return logprob;
}

/** Parameter settings
 *
 * See PARAM_* macros in MeanField.h for options
 */
void MeanField::setParameters (int param, void* value)
{

  switch (param)
    {
    case PARAM_DAMPER:
      {
        FloatType damp = *((FloatType*)value);

        if (damp>1 || damp<0)
          std::cerr << "MeanField::setParameters: PARAM_DAMPER must be " <<
            "between zero and one. Ignoring.\n";
        else
          m_beliefDamper = *((FloatType*)value);

        break;
      }
    case PARAM_BELDFUN:
      m_beliefDiffFun = *((DiffType*)value);
      break;

    case PARAM_BELDTOL:
      m_beliefDiffTol = *((FloatType*)value);
      break;

    default:
      std::cerr << "MeanField::setParameters: Unknown parameter index " <<
        param << '\n';
    }

}

 /** Add an edge to the graph.
 *
 * Note: Required by MRF, but NOT SUPPORTED
 */
/*
void MeanField::setNeighbors (int pix1, int pix2, CostVal weight)
{
  std::cerr << "MeanField::setNeighbors not supported\n";
  assert(false);
}
*/

/** Variational free energy  E_Q[U+V] - H(Q)
 *
 * Expected energy under the variational distribution minus the entropy of the
 * variational distribution
 */
FloatType MeanField::freeEnergy()
{
  register FloatType energy = 0;

  static int valll = 0;
  valll++;

  int pixel,pixel1;

  FloatType *belief, *belief1;

  if (m_grid_graph) {

	  for (int r=0 ; r<m_height ; r++)
		for (int c=0 ; c<m_width ; c++)
		  {
			pixel = pixelIndex(r,c);

			// Data terms
			belief = m_expNodeBeliefs[pixel];

			for (int b=0 ; b<m_nLabels ; b++ )
			  energy += belief[b] * getDataCost(pixel,b);

			// Smoothness terms

			if (r<m_height-1)
			  { // Down
				pixel1 = pixelIndex(r+1,c);

				belief = m_nodeBeliefs[pixel];
				belief1 = m_nodeBeliefs[pixel1];

				for (int b=0 ; b<m_nLabels ; b++ )
				  for (int b1=0 ; b1<m_nLabels ; b1++)
					energy += exp(belief[b]+belief1[b1]) *
					  getSmoothnessCost(pixel,pixel1,b,b1);

			  }

			if (c<m_width-1)
			  { // Right
				pixel1 = pixelIndex(r,c+1);

				belief = m_nodeBeliefs[pixel];
				belief1 = m_nodeBeliefs[pixel1];

				for (int b=0 ; b<m_nLabels ; b++ )
				  for (int b1=0 ; b1<m_nLabels ; b1++)
					energy += exp(belief[b]+belief1[b1]) *
					  getSmoothnessCost(pixel,pixel1,b,b1);
			  }
		  }

  } else { //general graphs

	  FloatType energySm = 0;


	  m_log->precision(2);
	  m_log->setf(std::ios::fixed,std::ios::floatfield);


	  for (int pixel=0 ; pixel<m_nPixels ; pixel++)
	  {
		// Data terms
		belief = m_expNodeBeliefs[pixel];

		float belb=0.0, gDC=0.0;

		for (int b=0 ; b<m_nLabels ; b++ ) {

//			if (valll==4)
//			*m_log <<" belief " << b << "  "<< belief[b] <<"  "<< getDataCost(pixel, b) <<endl;
			 belb = belief[b];
			 gDC = getDataCost(pixel,b);

			 // the condition below is to avoid mul(0, inf) which gives NaN.
			 if (belb==0.0)
				 energy += 0.0;
			 else
				 energy += belb * gDC;
		}

		// Smoothness terms

		// need to get neighbors

    	for (int p=0; p<m_nNeighbors[pixel]; p++) {

    		int pixel1=findNeighborWithIdx(pixel, p);
    		if (pixel1==-1) {
    			printf("Error\n");
    			exit(1);
    		}

			belief = m_nodeBeliefs[pixel];
			belief1 = m_nodeBeliefs[pixel1];

			for (int b=0 ; b<m_nLabels ; b++ )
			  for (int b1=0 ; b1<m_nLabels ; b1++)
				energySm += exp(belief[b]+belief1[b1]) *
				  getSmoothnessCost(pixel,pixel1,b,b1);

    	}
    	energy += energySm/2;

	  }
  }

  return energy - getTotalNodeEntropy();
}


/** Entropy of the variational distribution
*/
FloatType MeanField::getTotalNodeEntropy()
{
  register FloatType entropy = 0;

  FloatType** beliefs = m_nodeBeliefs;
  FloatType** ebeliefs = m_expNodeBeliefs;

  FloatType *bel,*ebel;

  for (int i=0 ; i<m_nPixels ; i++, beliefs++, ebeliefs++)
    {
      bel = *beliefs;
      ebel = *ebeliefs;

      for (int b=0 ; b<m_nLabels ; b++, bel++, ebel++) {
    	  if (*ebel == 0.0 || *bel == 0.0)
    		  entropy += 0.0;
    	  else
    		  entropy += (*ebel) * (*bel);
      }
    }

  return -entropy;
}

/** Entropy of an individual node's variational distribution
*/
FloatType MeanField::getNodeEntropy(const int pixel)
{
  register FloatType entropy = 0;

  FloatType* belief = m_nodeBeliefs[pixel];
  FloatType* ebelief = m_expNodeBeliefs[pixel];

  for (int b=0 ; b<m_nLabels ; b++, belief++, ebelief++)
        entropy += (*ebelief) * (*belief);


  return -entropy;
}

