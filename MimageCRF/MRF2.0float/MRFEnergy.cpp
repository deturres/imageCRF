#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "MRFEnergy.h"
#include "mrf.h"

// Data-value accessor macro
#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]

// Homogeneous smoothness-value accessor macro
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]

MRFEnergy::MRFEnergy(int width, int height, int nLabels, EnergyFunction *eng) :
  MRF(width,height,nLabels,eng)
{
  initMRFEnergy();
}


MRFEnergy::MRFEnergy(int nPixels, int nLabels,EnergyFunction *eng) :
  MRF(nPixels,nLabels,eng)
{
  initMRFEnergy();

  //chet - for non-grid-graphs
  m_neighbors = (LinkedBlockList *) new LinkedBlockList[nPixels];
  m_varWeights = false;
  m_nNeighbors = new int[m_nPixels];
}

MRFEnergy::~MRFEnergy ()
{
  delete[] m_labels;

  if (m_freeV)
    delete[] m_V;

  if (!m_grid_graph) {
	  delete[] m_nNeighbors;
	  delete[] m_neighbors;
  }


}

/* Common initialization routine */
void MRFEnergy::initMRFEnergy ()
{
  m_log = &(std::cout);
  m_freeV = false;
  m_labels = new Label[m_nPixels];
  m_imWriter = 0;
}



/****************************************
 * Total Energy Functions
 ****************************************/

/** Return the total energy for the data terms for the current assignment
 *
 *  Note: Assumes there is something meaningful in the answer array.
 */
MRF::EnergyVal MRFEnergy::dataEnergy ()
{

  EnergyVal eng = (EnergyVal) 0;


  if ( m_dataType == ARRAY)
    {
      for ( int i = 0; i < m_nPixels; i++ )
      {
        eng = eng + m_D(i,m_labels[i]);
      }
    }
  else
    {
      for ( int i = 0; i < m_nPixels; i++ )
        eng = eng + m_dataFn(i,m_labels[i]);
    }
  return(eng);

}

/** Return the total energy for the smoothness terms for the current assignment.
 *
 *  Note: Assumes there is something meaningful in the answer array.
 */
MRF::EnergyVal MRFEnergy::smoothnessEnergy ()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;

    int r,c,pix;

    if ( m_grid_graph )
	{

		if ( m_smoothType != FUNCTION  )
		  {
			// Horizontal edges
			for ( r = 0; r < m_height; r++ )
			  for ( c = 1; c < m_width; c++ )
				{
				  pix    = c+r*m_width;
				  weight = m_varWeights ? m_horizWeights[pix-1] :  1;
				  eng = eng + m_V(m_labels[pix], m_labels[pix-1]) * weight;
				}

			// Vertical edges
			for ( r = 1; r < m_height; r++ )
			  for ( c = 0; c < m_width; c++ )
				{
				  pix = c+r*m_width;
				  weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
				  eng = eng + m_V(m_labels[pix], m_labels[pix-m_width]) * weight;
				}
		  }
		else
		  {
			// Horizontal edges
			for ( r = 0; r < m_height; r++ )
			  for ( c = 1; c < m_width; c++ )
				{
				  pix = c+r*m_width;
				  eng = eng + m_smoothFn(pix, pix-1, m_labels[pix], m_labels[pix-1]);
				}

			// Vertical edges
			for ( r = 1; r < m_height; r++ )
			  for ( c = 0; c < m_width; c++ )
				{
				  pix = c+r*m_width;
				  eng = eng + m_smoothFn(pix, pix-m_width, m_labels[pix],
										 m_labels[pix-m_width]);
				}
		  }
	}
    else // general graph
	{
	    if ( m_smoothType != FUNCTION  )
		{
		    for (int i=0; i<m_nPixels; i++)
			{
		    	LinkedBlockList *m_NBTR = &m_neighbors[i];
		    	m_NBTR->setCursorFront();
		    	while(m_NBTR->hasNext()) {

		    		Neighbor *temp1 = (Neighbor *) m_NBTR->next();
		    	    int pix2  = temp1->to_node;

		    		eng = eng + m_V(m_labels[i],m_labels[pix2]);
		    	}
			}
		    eng = eng/2;
		}
	    else // smooth is a function
		{
		    for (int i=0; i<m_nPixels; i++)
			{
		    	LinkedBlockList *m_NBTR = &m_neighbors[i];
		    	m_NBTR->setCursorFront();
		    	while(m_NBTR->hasNext()) {
		    		Neighbor *temp1 = (Neighbor *) m_NBTR->next();
		    	    int pix2  = temp1->to_node;
		    		eng = eng + m_smoothFn(i, pix2, m_labels[i], m_labels[pix2]);
		    	}
			}
		    eng = eng/2;
		}
	}

	return(eng);

}

/*****************************************
 * Component costs of energy functions
 *****************************************/

/* Data cost from the energy function for a particular assignment
 */
MRF::CostVal MRFEnergy::getDataCost (const int pixel, const Label label)
{
  return (m_dataType == ARRAY ?
          m_D(pixel,label) :
          m_dataFn(pixel,label) );
}

/** Smoothness cost from the energy function for a particular assignment
 *
 * NOTE: Assumes a smoothness function handle handles symmetry.
 */
MRF::CostVal MRFEnergy::getSmoothnessCost(int pixel1, int pixel2,
                                                 Label label1, Label label2)
{

  // Order the pixel naming by value
  if (pixel1 > pixel2)
    {
      int tmpP = pixel1;
      pixel1 = pixel2;
      pixel2 = tmpP;

      Label tmpL = label1;
      label1 = label2;
      label2 = tmpL;
    }

  assert(pixel1<pixel2);

  if (m_smoothType == ARRAY)
    if (pixel1+1 == pixel2 && m_width>1) // Horizontal edge
      return
        (m_varWeights ? m_horizWeights[pixel1] : 1) * m_V(label1,label2) ;
    else
      return
        (m_varWeights ? m_vertWeights[pixel1] : 1) * m_V(label1,label2);
  else
    return m_smoothFn(pixel1,pixel2,label1,label2);



}

/****************************************
 * Aggregate cost functions
 ****************************************/

/** Copy energy terms for one pixel into a destination array
 */
void MRFEnergy::copyDataCost (const int pixel, CostVal* dst)
{
  if (m_dataType == ARRAY)
    {
      CostVal* pc = m_D + pixel*m_nLabels;

      for (int i=0 ; i<m_nLabels ; i++, pc++, dst++)
        *dst = - *pc;
    }
  else
    {
      for (int i=0 ; i<m_nLabels ; i++,dst++)
        *dst = - getDataCost(pixel, i);
    }
}


/* Add smoothness costs for a pair of neighboring pixels with one
 * fixed pixel label to a destination array.
 */
void MRFEnergy::addSmoothnessCost (const int fixedPixel,
                                            const int varPixel,
                                            const Label fixedLabel,
                                            CostVal* dst)
{

  for (int i=0 ; i<m_nLabels ; i++, dst++)
    {
      *dst -= getSmoothnessCost(fixedPixel,varPixel, fixedLabel, i);
    }

}

/** Add smoothness costs for a pair of neighboring pixels and write it
 *  to the destination array provided (must be of length m_nlabels^2).
 *
 * Note that pixel varies in the row and pixel2 varies in the column.
 */
void MRFEnergy::addSmoothnessCost(const int pixel1, const int pixel2,
                                    CostVal* dst)
{
  for (int i=0 ; i<m_nLabels ; i++)
    for (int j=0 ; j<m_nLabels ; j++, dst++)
      *dst -= getSmoothnessCost(pixel1,pixel2,i,j);
}

/** Copy sparse energy terms for one pixel into a destination array
 */
void MRFEnergy::copySparseDataCost (const int pixel, CostVal* dst, bool* sparsity)
{
  if (m_dataType == ARRAY)
    {
      CostVal* pc = m_D + pixel*m_nLabels;

      for (int i=0 ; i<m_nLabels ; i++, pc++, dst++,sparsity++)
        if (*sparsity)
          *dst = - *pc;
      // xxx SHOULD WE ZERO OTHERWISE?
    }
  else
    {
      for (int i=0 ; i<m_nLabels ; i++,dst++)
        *dst = - getDataCost(pixel, i);
    }
}

/* Add smoothness costs for an edge with one fixed pixel label to a destination
 * array
 */
void MRFEnergy::addSparseSmoothnessCost (const int fixedPixel,
                                                 const int varPixel,
                                                 const Label fixedLabel,
                                                 CostVal* dst, bool* sparsity)
{

  for (int i=0 ; i<m_nLabels ; i++, dst++, sparsity++)
    {
      if (*sparsity)
        *dst -= getSmoothnessCost(fixedPixel,varPixel, fixedLabel, i);
    }

}

/******************************************
 * Post-Optimization Label/Answer Functions
 *****************************************/

/** Returns a pointer to an array allowing for storage of the answer.
 *
 *  Note: This is somewhat unsafe, but required by MRF.
 */
MRF::Label* MRFEnergy::getAnswerPtr ()
{
  return m_labels;
}


/** Current answer at a pixel
 */
MRF::Label MRFEnergy::getLabel (int pixel)
{
  return m_labels[pixel];
}

/** Set the current answer for a pixel
 */
void MRFEnergy::setLabel (int pixel, Label label)
{
  m_labels[pixel] = label;
}

/** Use a function to write current answers to an image */
void MRFEnergy::setLabelImageWriter(LabelImageWriter* imWriter)
{
  m_imWriter = imWriter;
}

/** Clear current answers */
void MRFEnergy::clearAnswer()
{
  for (int i=0 ; i<m_nPixels ; i++)
    m_labels[i] = 0;
}

/****************************************
 * Cost Function Set Methods
 ****************************************/
#ifdef ENERGY_FRIEND_COSTS // This condition shouldn't be met probably
void MRFEnergy::setData(DataCost* dcost)
{

  if (dcost==0)
    {
      std::cerr << "WARNING: setData received null DataCost. Ignoring.";
      return;
    }

  if (dcost->m_type != m_dataType)
    {
      std::cerr << "WARNING: Cannot change data type. Ignoring setData.";
      return;
    }

  if (m_dataType == ARRAY)
    setData(dcost->m_costArray);
  else
    setData(dcost->m_costFn);

}

void MRFEnergy::setSmoothness(SmoothnessCost* scost)
{

  if (scost==0)
    {
      std::cerr << "WARNING: setSmoothness received null SmoothnessCost. " <<
        "Ignoring.";
      return;
    }

  if (scost->m_type != m_smoothType)
    {
      std::cerr << "WARNING: Cannot change smoothness type. " <<
        "Ignoring setSmoothness.";
      return;
    }

  if ( m_smoothType == FUNCTION )
    setSmoothness(scost->m_costFn);
  else
    {
      if ( m_smoothType == ARRAY )
        {
          checkArray(scostost->m_V);
          setSmoothness(scost->m_V);
        }
      else
        {
          int smoothExp = scost->m_smoothExp;
          if ( smoothExp != 1 && smoothExp != 2)
            {
              std::cerr << "WARNING: Wrong exponent in setSmoothness. " <<
                "Ignoring.";
              return;
            }
          setSmoothness(smoothExp,scost->m_smoothMax,scost->m_lambda);
        }

      if (scost->m_varWeights )
        {
          if (!m_grid_graph)
            std::cerr << "WARNING: Edge multiplier cannot be used with " <<
              "non-grid graphs. Ignoring.";
          else
            setCues(scost->m_hWeights,scostost->m_vWeights);
        }

    }

}
#endif

void MRFEnergy::setData(DataCostFn dcost)
{
    m_dataFn = dcost;
    m_dataType = FUNCTION;
}

void MRFEnergy::setData(CostVal* data)
{
  m_D = data;
  m_dataType = ARRAY;
}

void MRFEnergy::setSmoothness(SmoothCostGeneralFn cost)
{
  m_smoothFn = cost;
  m_smoothType = FUNCTION;
}

void MRFEnergy::setSmoothness(CostVal* V)
{
  if (m_freeV)
    {
      delete[] m_V;
      m_freeV = false;
    }

  m_V = V;
  m_smoothType = ARRAY;

}

/** Set a smoothness function of the form
 *    V(l1,l2) = lambda * min ( |l1-l2|^m_smoothExp, m_smoothMax )
 */
void MRFEnergy::setSmoothness(int smoothExp, CostVal smoothMax, CostVal lambda)
{

  m_smoothType = ARRAY;
  CostVal cost;

  if (!m_freeV) // space for the smoothness cost hasn't already been allocated
    {
      m_freeV = true;

      m_V = (CostVal *) new CostVal[m_nLabels*m_nLabels];
    }

  if (!m_V)
    {
      std::cerr << "Not enough memory!\n";
      exit(1);
    }


  for (int i=0; i<m_nLabels; i++)
    for (int j=i; j<m_nLabels; j++)
      {
        // L1 or L2 cost
        cost = (CostVal)((smoothExp == 1) ? j - i : (j - i)*(j - i));

        // Clip at smoothMax
        if (cost > smoothMax)
          cost = smoothMax;

        // Set at symmetric locations
        m_V[i*m_nLabels + j] = m_V[j*m_nLabels + i] = cost*lambda;
      }
}


/** Set spatially varying multipliers of the homogenous smoothness costs
 */
void MRFEnergy::setCues(CostVal* hCue, CostVal* vCue)
{
	m_horizWeights = hCue;
	m_vertWeights  = vCue;
    m_varWeights = true;
}

/** Access spatially varying multipliers of the homogenous smoothness costs
 *  for horizontal edges
 */
MRF::CostVal MRFEnergy::getHorizWeight(const int r, const int c)
{
  return m_varWeights ? m_horizWeights[r*m_width + c] :  1;
}

/** Access spatially varying multipliers of the homogenous smoothness costs
 *  for vertical edges
 */
MRF::CostVal MRFEnergy::getVertWeight(int r, int c)
{
  return  m_varWeights ? m_vertWeights[r*m_width + c] :  1;
}


/** Set an output stream for logging */
void MRFEnergy::setLog(std::ostream *stream) { m_log = stream; }


////chet functions added after this line

void MRFEnergy::setNeighbors(int pixel1, int pixel2, CostVal weight)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0 && pixel2 < m_nPixels && pixel2 >= 0);
    assert(m_grid_graph == 0);

    Neighbor *temp1 = (Neighbor *) new Neighbor;
    Neighbor *temp2 = (Neighbor *) new Neighbor;


    temp1->weight  = weight;
    temp1->to_node = pixel2;

    temp2->weight  = weight;
    temp2->to_node = pixel1;

//chet changed ttwo function calls below to use one with two parameters
// so that i can set m_neighborFlag to indicate type of pointer
// which is useful while deleting the pointer
    m_neighbors[pixel1].addFront(temp1, 1);
    m_neighbors[pixel2].addFront(temp2, 1);

}

int MRFEnergy::traverseNeighbors(int pixel1)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[pixel1];

	int cnter = 0;
	m_NBTR->setCursorFront();

	while(m_NBTR->hasNext()) {
		cnter++;
		m_NBTR->next();
	}

	return cnter;

//	printf("Pixel %d has neighbor count =  %d \n", pixel1, cnter);

}




int MRFEnergy::findIdxOfNeighbor(int mainpix, int pixel1)
{

    assert(pixel1 < m_nPixels && pixel1 >= 0);
    assert(mainpix < m_nPixels && mainpix >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[mainpix];

	m_NBTR->setCursorFront();

	int idx=-1;

	int cnter=0;
	while(m_NBTR->hasNext()) {
		Neighbor* nb = (Neighbor *)m_NBTR->next();
		if (nb->to_node == pixel1) {
			idx = cnter;
			break;
		}
		cnter++;
	}

	return idx;
}

int MRFEnergy::findNeighborWithIdx(int mainpix, int idx)
{

    assert(mainpix < m_nPixels && mainpix >= 0);
    assert(m_grid_graph == 0);

	LinkedBlockList *m_NBTR = &m_neighbors[mainpix];

	int pixel=-1;

	int cnter = 0;
	m_NBTR->setCursorFront();

	while(m_NBTR->hasNext()) {
		if (idx==cnter) {
			Neighbor* nb = (Neighbor *) m_NBTR->next();
			pixel = nb->to_node;
			break;
		} else {
			cnter++;
			m_NBTR->next();
		}
	}

	return pixel;
}



void MRFEnergy::fillAnswersK() {
	for (int i=0; i<m_nPixels; i++) {
		m_labels[i] = m_nLabels;
	}
}


//if (m_grid_graph)

