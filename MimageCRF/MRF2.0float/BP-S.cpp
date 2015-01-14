#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <new>
#include "BP-S.h"

#define private public
#include "typeTruncatedQuadratic2D.h"
#undef private


#define m_D(pix,l)  m_D[(pix)*m_nLabels+(l)]
#define m_V(l1,l2)  m_V[(l1)*m_nLabels+(l2)]


#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define TRUNCATE_MIN(a,b) { if ((a) > (b)) (a) = (b); }
#define TRUNCATE_MAX(a,b) { if ((a) < (b)) (a) = (b); }
#define TRUNCATE TRUNCATE_MIN

/////////////////////////////////////////////////////////////////////////////
//                  Operations on vectors (arrays of size K)               //
/////////////////////////////////////////////////////////////////////////////

inline void CopyVector(BPS::REAL* to, MRF::CostVal* from, int K)
{
    BPS::REAL* to_finish = to + K;
    do
	{
	    *to ++ = *from ++;
	} while (to < to_finish);
}

inline void AddVector(BPS::REAL* to, BPS::REAL* from, int K)
{
    BPS::REAL* to_finish = to + K;
    do
	{
	    *to ++ += *from ++;
	} while (to < to_finish);
}

inline BPS::REAL SubtractMin(BPS::REAL *D, int K)
{
    int k;
    BPS::REAL delta;

    delta = D[0];
    for (k=1; k<K; k++) TRUNCATE(delta, D[k]);
    for (k=0; k<K; k++) D[k] -= delta;

    return delta;
}

// Functions UpdateMessageTYPE (see the paper for details):
//
// - Set Di[ki] := gamma*Di_hat[ki] - M[ki]
// - Set M[kj] := min_{ki} (Di[ki] + V[ki,kj])
// - Normalize message:
//        delta := min_{kj} M[kj]
//        M[kj] := M[kj] - delta
//        return delta
//
// If dir = 1, then the meaning of i and j is swapped.

///////////////////////////////////////////
//                  L1                   //
///////////////////////////////////////////

inline BPS::REAL UpdateMessageL1(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal smoothMax)
{
    int k;
    BPS::REAL delta;

    delta = M[0] = gamma*Di_hat[0] - M[0];
    for (k=1; k<K; k++)
	{
	    M[k] = gamma*Di_hat[k] - M[k];
	    TRUNCATE(delta, M[k]);
	    TRUNCATE(M[k], M[k-1] + lambda);
	}

    M[--k] -= delta;
    TRUNCATE(M[k], lambda*smoothMax);
    for (k--; k>=0; k--)
	{
	    M[k] -= delta;
	    TRUNCATE(M[k], M[k+1] + lambda);
	    TRUNCATE(M[k], lambda*smoothMax);
	}

    return delta;
}

////////////////////////////////////////
//               L2                   //
////////////////////////////////////////

inline BPS::REAL UpdateMessageL2(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal smoothMax, void *buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int* parabolas = (int*) ((char*)buf + K*sizeof(BPS::REAL));
    int* intersections = parabolas + K;
    TypeTruncatedQuadratic2D::REAL* Di_tmp = (TypeTruncatedQuadratic2D::REAL*) (intersections + K + 1);
    TypeTruncatedQuadratic2D::REAL* M_tmp = Di_tmp + K;
    TypeTruncatedQuadratic2D::Edge* tmp = NULL;

    int k;
    BPS::REAL delta;

    assert(lambda >= 0);

    Di[0] = gamma*Di_hat[0] - M[0];
    delta = Di[0];
    for (k=1; k<K; k++)
	{
	    Di[k] = gamma*Di_hat[k] - M[k];
	    TRUNCATE(delta, Di[k]);
	}

    if (lambda == 0)
	{
	    for (k=0; k<K; k++) M[k] = 0;
	    return delta;
	}

    for (k=0; k<K; k++) Di_tmp[k] = Di[k];
    tmp->DistanceTransformL2(K, 1, lambda, Di_tmp, M_tmp, parabolas, intersections);
    for (k=0; k<K; k++) M[k] = (BPS::REAL) M_tmp[k];

    for (k=0; k<K; k++)
	{
	    M[k] -= delta;
	    TRUNCATE(M[k], lambda*smoothMax);
	}

    return delta;
}


//////////////////////////////////////////////////
//                FIXED_MATRIX                  //
//////////////////////////////////////////////////

inline BPS::REAL UpdateMessageFIXED_MATRIX(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal* V, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = gamma*Di_hat[ki] - M[ki];
	}

    for (kj=0; kj<K; kj++)
	{
	    M[kj] = Di[0] + lambda*V[0];
	    V ++;
	    for (ki=1; ki<K; ki++)
		{
		    TRUNCATE(M[kj], Di[ki] + lambda*V[0]);
		    V ++;
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}

/////////////////////////////////////////////
//                GENERAL                  //
/////////////////////////////////////////////

inline BPS::REAL UpdateMessageGENERAL(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, int dir, MRF::CostVal* V, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = (gamma*Di_hat[ki] - M[ki]);
	}

    if (dir == 0)
	{
	    for (kj=0; kj<K; kj++)
		{
		    M[kj] = Di[0] + V[0];
		    V ++;
		    for (ki=1; ki<K; ki++)
			{
			    TRUNCATE(M[kj], Di[ki] + V[0]);
			    V ++;
			}
		}
	}
    else
	{
	    for (kj=0; kj<K; kj++)
		{
		    M[kj] = Di[0] + V[0];
		    V += K;
		    for (ki=1; ki<K; ki++)
			{
			    TRUNCATE(M[kj], Di[ki] + V[0]);
			    V += K;
			}
		    V -= K*K - 1;
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}

inline BPS::REAL UpdateMessageGENERAL(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, BPS::SmoothCostGeneralFn fn, int i, int j, void* buf)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = (gamma*Di_hat[ki] - M[ki]);
	}

    for (kj=0; kj<K; kj++)
	{
	    M[kj] = Di[0] + fn(i, j, 0, kj);
	    for (ki=1; ki<K; ki++)
		{
		    delta = Di[ki] + fn(i, j, ki, kj);
		    TRUNCATE(M[kj], delta);
		}
	}

    delta = M[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, M[kj]);
    for (kj=0; kj<K; kj++) M[kj] -= delta;

    return delta;
}






BPS::BPS(int width, int height, int nLabels,EnergyFunction *eng):MRF(width,height,nLabels,eng)
{
    Allocate();
}
BPS::BPS(int nPixels, int nLabels,EnergyFunction *eng):MRF(nPixels,nLabels,eng)
{
    Allocate();

    //chet - for non-grid-graphs
    m_neighbors = (LinkedBlockList *) new LinkedBlockList[nPixels];
    m_varWeights = false;

}

BPS::~BPS()
{
    delete[] m_answer;
    if ( m_needToFreeD ) delete [] m_D;
    if ( m_needToFreeV ) delete [] m_V;
    if ( m_messages ) delete [] m_messages;
    if ( m_DBinary ) delete [] m_DBinary;
    if ( m_horzWeightsBinary ) delete [] m_horzWeightsBinary;
    if ( m_vertWeightsBinary ) delete [] m_vertWeightsBinary;

    if (!m_grid_graph) {
        // need to delete m_messagesNG
        // need to delete m_neighbors
    	for (int i=0; i<m_nPixels; i++) {
    		delete[] m_messagesNG[i];
    	}
    	delete[] m_messagesNG;

    	delete[] m_neighbors;
    	delete[] m_nNeighbors;
    }


}


void BPS::Allocate()
{
    m_type = NONE;
    m_needToFreeV = false;
    m_needToFreeD = false;

    m_D = NULL;
    m_V = NULL;
    m_horzWeights = NULL;
    m_vertWeights = NULL;
    m_horzWeightsBinary = NULL;
    m_vertWeightsBinary = NULL;

    m_DBinary = NULL;
    m_messages = NULL;
    m_messageArraySizeInBytes = 0;

    m_answer = new Label[m_nPixels];

    m_messagesNG = NULL;
}

void BPS::clearAnswer()
{
    memset(m_answer, 0, m_nPixels*sizeof(Label));

    if (m_grid_graph) {
		if (m_messages)
		{
			memset(m_messages, 0, m_messageArraySizeInBytes);
		}
    } //else we have already set each to 0 while creating the m_messageNG structure.
}


MRF::EnergyVal BPS::smoothnessEnergy()
{
    EnergyVal eng = (EnergyVal) 0;
    EnergyVal weight;
    int x,y,pix;

    if ( m_grid_graph )
	{
	    if ( m_smoothType != FUNCTION  )
		{
		    for ( y = 0; y < m_height; y++ )
			for ( x = 1; x < m_width; x++ )
			    {
				pix    = x+y*m_width;
				weight = m_varWeights ? m_horzWeights[pix-1] :  1;
				eng = eng + m_V(m_answer[pix],m_answer[pix-1])*weight;
			    }

		    for ( y = 1; y < m_height; y++ )
			for ( x = 0; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				weight = m_varWeights ? m_vertWeights[pix-m_width] :  1;
				eng = eng + m_V(m_answer[pix],m_answer[pix-m_width])*weight;
			    }
		}
	    else
		{
		    for ( y = 0; y < m_height; y++ )
			for ( x = 1; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				eng = eng + m_smoothFn(pix,pix-1,m_answer[pix],m_answer[pix-1]);
			    }

		    for ( y = 1; y < m_height; y++ )
			for ( x = 0; x < m_width; x++ )
			    {
				pix = x+y*m_width;
				eng = eng + m_smoothFn(pix,pix-m_width,m_answer[pix],m_answer[pix-m_width]);
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

		    		eng = eng + m_V(m_answer[i],m_answer[pix2]);
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
		    		eng = eng + m_smoothFn(i, pix2, m_answer[i], m_answer[pix2]);
		    	}
			}
		    eng = eng/2;
		}
	}

    return(eng);
}


// this is local cost
MRF::EnergyVal BPS::dataEnergy()
{
  //  EnergyVal eng = (EnergyVal) 0;
  double eng = 0;

  if ( m_dataType == ARRAY)
	{
	  for ( int i = 0; i < m_nPixels; i++ )
	    eng = eng + m_D(i,m_answer[i]);
	}
  else
	{
    for ( int i = 0; i < m_nPixels; i++ )
      eng = eng + m_dataFn(i,m_answer[i]);
	}

  return((EnergyVal)eng);

}


void BPS::setData(DataCostFn dcost)
{
    int i, k;

    m_dataFn = dcost;
    CostVal* ptr;
    m_D = new CostVal[m_nPixels*m_nLabels];

    for (ptr=m_D, i=0; i<m_nPixels; i++)
	for (k=0; k<m_nLabels; k++, ptr++)
	    {
		*ptr = m_dataFn(i,k);
	    }
    m_needToFreeD = true;
}

void BPS::setData(CostVal* data)
{
    m_D = data;
    m_needToFreeD = false;
}


void BPS::setSmoothness(SmoothCostGeneralFn cost)
{
    assert(m_horzWeights == NULL && m_vertWeights == NULL && m_V == NULL);

    int x, y, i, ki, kj;
    CostVal* ptr;

    m_smoothFn = cost;
    m_type = GENERAL;

    if (!m_allocateArrayForSmoothnessCostFn) return;

    return; //chet forcing the use of function instead of array
    // the code below cannot be used for general graphs.
    // anyway everything below is a efficiency hack

    // try to cache all the function values in an array for efficiency
    m_V = new(std::nothrow) CostVal[2*m_nPixels*m_nLabels*m_nLabels];
    if (!m_V) {
	fprintf(stderr, "not caching smoothness cost values (not enough memory)\n");
	return; // if not enough space, just call the function directly
    }

    m_needToFreeV = true;

    for (ptr=m_V,i=0,y=0; y<m_height; y++)
	for (x=0; x<m_width; x++, i++)
	    {
		if (x < m_width-1)
		    {
			for (kj=0; kj<m_nLabels; kj++)
			    for (ki=0; ki<m_nLabels; ki++)
				{
				    *ptr++ = cost(i,i+1,ki,kj);
				}
		    }
		else ptr += m_nLabels*m_nLabels;

		if (y < m_height-1)
		    {
			for (kj=0; kj<m_nLabels; kj++)
			    for (ki=0; ki<m_nLabels; ki++)
				{
				    *ptr++ = cost(i,i+m_width,ki,kj);
				}
		    }
		else ptr += m_nLabels*m_nLabels;
	    }
}
void BPS::setSmoothness(CostVal* V)
{
    m_type = FIXED_MATRIX;
    m_V = V;
}


void BPS::setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda)
{
    assert(smoothExp == 1 || smoothExp == 2);
    assert(lambda >= 0);

    m_type = (smoothExp == 1) ? L1 : L2;

    int ki, kj;
    CostVal cost;

    m_needToFreeV = true;

    m_V = new CostVal[m_nLabels*m_nLabels];

    for (ki=0; ki<m_nLabels; ki++)
	for (kj=ki; kj<m_nLabels; kj++)
	    {
		cost = (CostVal) ((smoothExp == 1) ? kj - ki : (kj - ki)*(kj - ki));
		if (cost > smoothMax) cost = smoothMax;
		m_V[ki*m_nLabels + kj] = m_V[kj*m_nLabels + ki] = cost*lambda;
	    }

    m_smoothMax = smoothMax;
    m_lambda = lambda;
}


void BPS::setCues(CostVal* hCue, CostVal* vCue)
{
    m_horzWeights = hCue;
    m_vertWeights  = vCue;
}


void BPS::initializeAlg()
{
    assert(m_type != NONE);

    if (m_grid_graph) {

		int i;

		// determine type
		if (m_type == L1 && m_nLabels == 2)
		{
			m_type = BINARY;
		}

		// allocate messages
		int messageNum = (m_type == BINARY) ? 4*m_nPixels : 4*m_nPixels*m_nLabels;
		m_messageArraySizeInBytes = messageNum*sizeof(REAL);
		m_messages = new REAL[messageNum];
		memset(m_messages, 0, messageNum*sizeof(REAL));

		if (m_type == BINARY)
		{
			assert(m_DBinary == NULL && m_horzWeightsBinary == NULL && m_horzWeightsBinary == NULL);
			m_DBinary = new CostVal[m_nPixels];
			m_horzWeightsBinary = new CostVal[m_nPixels];
			m_vertWeightsBinary = new CostVal[m_nPixels];

			if ( m_dataType == ARRAY)
			{
				for (i=0; i<m_nPixels; i++)
				{
					m_DBinary[i] = m_D[2*i+1] - m_D[2*i];
				}
			}
			else
			{
				for (i=0; i<m_nPixels; i++)
				{
					m_DBinary[i] = m_dataFn(i,1) - m_dataFn(i,0);
				}
			}

			assert(m_V[0] == 0 && m_V[1] == m_V[2] && m_V[3] == 0);
			for (i=0; i<m_nPixels; i++)
			{
				m_horzWeightsBinary[i] = (m_varWeights) ? m_V[1]*m_horzWeights[i] : m_V[1];
				m_vertWeightsBinary[i] = (m_varWeights) ? m_V[1]*m_vertWeights[i] : m_V[1];
			}
		}

    } else { //general graphs

        // allocate messages
        // this is message number per neighbor
        int messageNumPerNB = m_nPixels*m_nLabels;
        // true number of messages will be this multiplied with the number of neighbors per pixel-label
        // but in general graph this number of neighbors could be variable

        m_messagesNG = new REAL*[m_nPixels];
        m_nNeighbors = new int[m_nPixels];

        // need to traverse the neighbor list to get neighbors per pixel

        int totNB=0;

        m_maxnNB=0;

    	// traverse Neighbors list
    	for (int m=0; m < m_nPixels; m++) {
    		int numNB=traverseNeighbors((int)m);

    		m_nNeighbors[m] = numNB;

    		if (m_maxnNB<numNB) {
    			m_maxnNB=numNB;
    		}

    		m_messagesNG[m] = new REAL[numNB*m_nLabels];
    		memset(m_messagesNG[m], 0, numNB*m_nLabels*sizeof(REAL));

    		totNB +=numNB;
    	}

        m_messageArraySizeInBytes = totNB*sizeof(REAL);
    }

}

void BPS::optimizeAlg(int nIterations)
{
    assert(m_type != NONE);

    if (m_grid_graph)
	{
	    switch (m_type)
		{
		case L1:            optimize_GRID_L1(nIterations);           break;
		case L2:            optimize_GRID_L2(nIterations);           break;
		case FIXED_MATRIX:  optimize_GRID_FIXED_MATRIX(nIterations); break;
		case GENERAL:       optimize_GRID_GENERAL(nIterations);      break;
		case BINARY:        optimize_GRID_BINARY(nIterations);       break;
		default: assert(0); exit(1);
		}
	}
    else {  // general graph

	    switch (m_type)
		{
		case FIXED_MATRIX:
		case GENERAL:       optimize_NONGRID_FIXED_MATRIX_GENERAL(nIterations);      break;
		default: assert(0); exit(1);
		}


    }

    //	printf("lower bound = %f\n", m_lowerBound);

    ////////////////////////////////////////////////
    //          computing solution                //
    ////////////////////////////////////////////////

    // fill answers with m_nLabel so that we can know which nodes are labeled
     // and which are not while computeSolution fills the nodes.
     fillAnswersK();

     if (m_grid_graph)
 	{
 	    switch (m_type)
 		{
 		case L1:
 		case L2:
 		case FIXED_MATRIX:  computeSolution_GRID_L1_L2_FIXED_MATRIX(); break;
 		case GENERAL:       computeSolution_GRID_GENERAL();      break;
 		case BINARY:        computeSolution_GRID_BINARY();       break;
 		default: assert(0); exit(1);
 		}
 	}
     else {  // general graph

 	    switch (m_type)
 		{
 		case FIXED_MATRIX:
 		case GENERAL:       computeSolution_NONGRID_FIXED_MATRIX_GENERAL();      break;
 		default: assert(0); exit(1);
 		}
     }

/*

    if (m_type != BINARY)
	{
	    int x, y, n, K = m_nLabels;
	    CostVal* D_ptr;
	    REAL* M_ptr;
	    REAL* Di;
	    REAL delta;
	    int ki, kj;

	    Di = new REAL[K];

	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);

			if (m_type == GENERAL)
			    {
				if (m_V)
				    {
					CostVal* ptr = m_V + 2*(x+y*m_width-1)*K*K;
					if (x > 0)
					    {
						kj = m_answer[n-1];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += ptr[kj + ki*K];
						    }
					    }
					ptr -= (2*m_width-3)*K*K;
					if (y > 0)
					    {
						kj = m_answer[n-m_width];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += ptr[kj + ki*K];
						    }
					    }
				    }
				else
				    {
					if (x > 0)
					    {
						kj = m_answer[n-1];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += m_smoothFn(n, n-1, ki, kj);
						    }
					    }
					if (y > 0)
					    {
						kj = m_answer[n-m_width];
						for (ki=0; ki<K; ki++)
						    {
							Di[ki] += m_smoothFn(n, n-m_width, ki, kj);
						    }
					    }
				    }
			    }
			else // m_type == L1, L2 or FIXED_MATRIX
			    {
				if (x > 0)
				    {
					kj = m_answer[n-1];
					CostVal lambda = (m_varWeights) ? m_horzWeights[n-1] : 1;
					for (ki=0; ki<K; ki++)
					    {
						Di[ki] += lambda*m_V[kj*K + ki];
					    }
				    }
				if (y > 0)
				    {
					kj = m_answer[n-m_width];
					CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
					for (ki=0; ki<K; ki++)
					    {
						Di[ki] += lambda*m_V[kj*K + ki];
					    }
				    }
			    }

			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			// compute min
			delta = Di[0];
			m_answer[n] = 0;
			for (ki=1; ki<K; ki++)
			    {
				if (delta > Di[ki])
				    {
					delta = Di[ki];
					m_answer[n] = ki;
				    }
			    }
		    }

	    delete [] Di;
	}
    else // m_type == BINARY
	{
	    int x, y, n;
	    REAL* M_ptr;
	    REAL Di;

	    n = 0;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, M_ptr+=2, n++)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += (m_answer[n-1] == 0)       ? m_horzWeightsBinary[n-1]       : -m_horzWeightsBinary[n-1];
			if (y > 0) Di += (m_answer[n-m_width] == 0) ? m_vertWeightsBinary[n-m_width] : -m_vertWeightsBinary[n-m_width];

			if (x < m_width-1)  Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			// compute min
			m_answer[n] = (Di >= 0) ? 0 : 1;
		    }
	}

	*/

}

void BPS::optimize_GRID_L1(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;

    Di = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n] : m_lambda;
				UpdateMessageL1(M_ptr, Di, K, 1, lambda, m_smoothMax);
			    }
			if (y < m_height-1)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n] : m_lambda;
				UpdateMessageL1(M_ptr+K, Di, K, 1, lambda, m_smoothMax);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n-1] : m_lambda;
				UpdateMessageL1(M_ptr-2*K, Di, K, 1, lambda, m_smoothMax);
			    }
			if (y > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n-m_width] : m_lambda;
				UpdateMessageL1(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_smoothMax);
			    }
		    }
	}

    delete [] Di;
}

void BPS::optimize_GRID_L2(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new char[2*K*sizeof(TypeTruncatedQuadratic2D::REAL) + (2*K+1)*sizeof(int) + K*sizeof(REAL)];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n] : m_lambda;
				UpdateMessageL2(M_ptr, Di, K, 1, lambda, m_smoothMax, buf);
			    }
			if (y < m_height-1)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n] : m_lambda;
				UpdateMessageL2(M_ptr+K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_horzWeights[n-1] : m_lambda;
				UpdateMessageL2(M_ptr-2*K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
			if (y > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_lambda*m_vertWeights[n-m_width] : m_lambda;
				UpdateMessageL2(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_smoothMax, buf);
			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}


void BPS::optimize_GRID_BINARY(int nIterations)
{
    int x, y, n;
    REAL* M_ptr;
    REAL Di;

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, M_ptr+=2, n++)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += M_ptr[-2]; // message (x-1,y)->(x,y)
			if (y > 0) Di += M_ptr[-2*m_width+1]; // message (x,y-1)->(x,y)
			if (x < m_width-1) Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			REAL DiScaled = Di * 1;
			if (x < m_width-1)
			    {
				Di = DiScaled - M_ptr[0];
				CostVal lambda = m_horzWeightsBinary[n];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[0] = lambda;
				else             M_ptr[0] = (Di < -lambda) ? -lambda : Di;
			    }
			if (y < m_height-1)
			    {
				Di = DiScaled - M_ptr[1];
				CostVal lambda = m_vertWeightsBinary[n];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[1] = lambda;
				else             M_ptr[1] = (Di < -lambda) ? -lambda : Di;
			    }
		    }

	    // backward pass
	    n --;
	    M_ptr -= 2;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, M_ptr-=2, n--)
		    {
			Di = m_DBinary[n];
			if (x > 0) Di += M_ptr[-2]; // message (x-1,y)->(x,y)
			if (y > 0) Di += M_ptr[-2*m_width+1]; // message (x,y-1)->(x,y)
			if (x < m_width-1) Di += M_ptr[0]; // message (x+1,y)->(x,y)
			if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

			REAL DiScaled = Di * 1;
			if (x > 0)
			    {
				Di = DiScaled - M_ptr[-2];
				CostVal lambda = m_horzWeightsBinary[n-1];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[-2] = lambda;
				else             M_ptr[-2] = (Di < -lambda) ? -lambda : Di;
			    }
			if (y > 0)
			    {
				Di = DiScaled - M_ptr[-2*m_width+1];
				CostVal lambda = m_vertWeightsBinary[n-m_width];
				if (lambda < 0) { Di = -Di; lambda = -lambda; }
				if (Di > lambda) M_ptr[-2*m_width+1] = lambda;
				else             M_ptr[-2*m_width+1] = (Di < -lambda) ? -lambda : Di;
			    }
		    }
	}
}

void BPS::optimize_GRID_FIXED_MATRIX(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1)
			    {
				CostVal lambda = (m_varWeights) ? m_horzWeights[n] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr, Di, K, 1, lambda, m_V, buf);
			    }
			if (y < m_height-1)
			    {
				CostVal lambda = (m_varWeights) ? m_vertWeights[n] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr+K, Di, K, 1, lambda, m_V, buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			SubtractMin(Di, K);

			if (x > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_horzWeights[n-1] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr-2*K, Di, K, 1, lambda, m_V, buf);
			    }
			if (y > 0)
			    {
				CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
				UpdateMessageFIXED_MATRIX(M_ptr-(2*m_width-1)*K, Di, K, 1, lambda, m_V, buf);
			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}

void BPS::optimize_GRID_GENERAL(int nIterations)
{
    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new REAL[K];

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    n = 0;
	    D_ptr = m_D;
	    M_ptr = m_messages;
	    CostVal* V_ptr = m_V;

	    for (y=0; y<m_height; y++)
		for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, V_ptr+=2*K*K, n++)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			if (x < m_width-1)
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr, Di, K, 1, /* forward dir*/ 0, V_ptr, buf);
				else     UpdateMessageGENERAL(M_ptr, Di, K, 1,   m_smoothFn, n, n+1,      buf);
			    }
			if (y < m_height-1)
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr+K, Di, K, 1, /* forward dir*/ 0, V_ptr+K*K, buf);
				else     UpdateMessageGENERAL(M_ptr+K, Di, K, 1,   m_smoothFn, n, n+m_width,    buf);
			    }
		    }

	    // backward pass
	    n --;
	    D_ptr -= K;
	    M_ptr -= 2*K;
	    V_ptr -= 2*K*K;

	    for (y=m_height-1; y>=0; y--)
		for (x=m_width-1; x>=0; x--, D_ptr-=K, M_ptr-=2*K, V_ptr-=2*K*K, n--)
		    {
			CopyVector(Di, D_ptr, K);
			if (x > 0) AddVector(Di, M_ptr-2*K, K); // message (x-1,y)->(x,y)
			if (y > 0) AddVector(Di, M_ptr-(2*m_width-1)*K, K); // message (x,y-1)->(x,y)
			if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
			if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

			// normalize Di, update lower bound
			SubtractMin(Di, K);

			if (x > 0)
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr-2*K, Di, K, 1, /* backward dir */ 1, V_ptr-2*K*K, buf);
				else     UpdateMessageGENERAL(M_ptr-2*K, Di, K, 1,   m_smoothFn, n, n-1,              buf);
			    }
			if (y > 0)
			    {
				if (m_V) UpdateMessageGENERAL(M_ptr-(2*m_width-1)*K, Di, K, 1, /* backward dir */ 1, V_ptr-(2*m_width-1)*K*K, buf);
				else     UpdateMessageGENERAL(M_ptr-(2*m_width-1)*K, Di, K, 1,   m_smoothFn, n, n-m_width,                    buf);

			    }
		    }
	}

    delete [] Di;
    delete [] (REAL *)buf;
}

//chet functions added after this line

void BPS::printMat(REAL* M, int K) {
	for (int i=0; i<K; i++) {
		printf("%5.2f  ", *M++);
	}
	printf("\n");
}

void BPS::setNeighbors(int pixel1, int pixel2, CostVal weight)
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

int BPS::traverseNeighbors(int pixel1)
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

void BPS::optimize_NONGRID_FIXED_MATRIX_GENERAL(int nIterations)
{
    int K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* N_ptr;
    REAL* Di;
    void* buf;

    Di = new REAL[K];
    buf = new REAL[K];

    CostVal lambda =  1;

    for ( ; nIterations > 0; nIterations --)
	{
	    // forward pass
	    D_ptr = m_D;

	    for (int m=0; m<m_nPixels; m++, D_ptr+=K) {
	    	CopyVector(Di, D_ptr, K);
//	    	printMat(Di, K);
	    	M_ptr = m_messagesNG[m];
	    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=K) {
	    		AddVector(Di, M_ptr, K);
//				printMat(M_ptr, K);
//				printMat(Di, K);
	    	}

	    	// update half of the messages in the forward pass
	    	M_ptr = m_messagesNG[m];
	    	for (int p=0; p<m_nNeighbors[m]/2; p++, M_ptr+=K) {
	    		// M_ptr should point to that incoming message
	    		// which edge has the outgoing message we are updating
	    		// in addition we need to pass the pointer to the outgoing message we are updating
	    		// hence need to find that out
	    		int q=findNeighborWithIdx(m, p);
	    		if (q==-1) {
	    			printf("Error\n");
	    			exit(1);
	    		}
	    		int idx = findIdxOfNeighbor(q, m);
	    		if (idx==-1) {
	    			printf("Error\n");
	    			exit(1);
	    		}
	    		N_ptr = m_messagesNG[q];
	    		N_ptr += idx*K;
	    		if ( m_smoothType != FUNCTION  ) { //assume ARRAY
	    			UpdateMessageFIXED_MATRIX_NG(M_ptr, Di, K, 1, lambda, m_V, buf, N_ptr);
	    		} else { // this is smoothfunction is function()
	    			UpdateMessageGENERAL_NG(M_ptr, Di, K, 1, lambda, m, q, buf, N_ptr);
	    		}
	    	}

	    }

	    // backward pass
	    D_ptr -= K;


	    for (int m=m_nPixels-1; m>=0; m--, D_ptr-=K) {

	    	CopyVector(Di, D_ptr, K);
//	    	printMat(Di, K);
	    	M_ptr = m_messagesNG[m];
	    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=K) {
	    		AddVector(Di, M_ptr, K);
//				printMat(M_ptr, K);
//				printMat(Di, K);
	    	}

	    	SubtractMin(Di, K);
	    	M_ptr = m_messagesNG[m];

	    	for (int p=m_nNeighbors[m]/2; p<m_nNeighbors[m]; p++, M_ptr+=K) {
	    		// M_ptr should point to that incoming message
	    		// which edge has the outgoing message we are updating
	    		// in addition we need to pass the pointer to the outgoing message we are updating
	    		// hence need to find that out
	    		int q=findNeighborWithIdx(m, p);
	    		if (q==-1) {
	    			printf("Error\n");
	    			exit(1);
	    		}
	    		int idx = findIdxOfNeighbor(q, m);
	    		if (idx==-1) {
	    			printf("Error\n");
	    			exit(1);
	    		}
	    		N_ptr = m_messagesNG[q];
	    		N_ptr += idx*K;

	    		if ( m_smoothType != FUNCTION  ) { //assume ARRAY
	    			UpdateMessageFIXED_MATRIX_NG(M_ptr, Di, K, 1, lambda, m_V, buf, N_ptr);
	    		} else { // this is smoothfunction is function()
	    			UpdateMessageGENERAL_NG(M_ptr, Di, K, 1, lambda, m, q, buf, N_ptr);
	    		}

	    	}

	    }

	} // iterations

    delete [] Di;
    delete [] (REAL *)buf;
}


int BPS::findIdxOfNeighbor(int mainpix, int pixel1)
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

int BPS::findNeighborWithIdx(int mainpix, int idx)
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

BPS::REAL BPS::UpdateMessageGENERAL_NG(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, int m, int q, void* buf, BPS::REAL* N)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = gamma*Di_hat[ki] - M[ki];
	}

    for (kj=0; kj<K; kj++)
	{
	    N[kj] = Di[0] + lambda*m_smoothFn(m, q, 0, kj);//[0];

	    for (ki=1; ki<K; ki++)
		{
		    TRUNCATE(N[kj], Di[ki] + lambda*m_smoothFn(m, q, ki, kj)); //V[0]);
		}
	}

    delta = N[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, N[kj]);
    for (kj=0; kj<K; kj++) N[kj] -= delta;

    return delta;
}


BPS::REAL BPS::UpdateMessageFIXED_MATRIX_NG(BPS::REAL* M, BPS::REAL* Di_hat, int K, BPS::REAL gamma, MRF::CostVal lambda, MRF::CostVal* V, void* buf, BPS::REAL* N)
{
    BPS::REAL* Di = (BPS::REAL*) buf;
    int ki, kj;
    BPS::REAL delta;

    for (ki=0; ki<K; ki++)
	{
	    Di[ki] = gamma*Di_hat[ki] - M[ki];
	}

    for (kj=0; kj<K; kj++)
	{
	    N[kj] = Di[0] + lambda*V[0];
	    V ++;
	    for (ki=1; ki<K; ki++)
		{
		    TRUNCATE(N[kj], Di[ki] + lambda*V[0]);
		    V ++;
		}
	}

    delta = N[0];
    for (kj=1; kj<K; kj++) TRUNCATE(delta, N[kj]);
    for (kj=0; kj<K; kj++) N[kj] -= delta;

    return delta;
}

void BPS::computeSolution_GRID_L1_L2_FIXED_MATRIX() {

    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    REAL delta;
    int ki, kj;

    Di = new REAL[K];

    n = 0;
    D_ptr = m_D;
    M_ptr = m_messages;

    for (y=0; y<m_height; y++)
	for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
	    {
		CopyVector(Di, D_ptr, K);

		if (x > 0)
		    {
			kj = m_answer[n-1];
			CostVal lambda = (m_varWeights) ? m_horzWeights[n-1] : 1;
			for (ki=0; ki<K; ki++)
			    {
				Di[ki] += lambda*m_V[kj*K + ki];
			    }
		    }
		if (y > 0)
		    {
			kj = m_answer[n-m_width];
			CostVal lambda = (m_varWeights) ? m_vertWeights[n-m_width] : 1;
			for (ki=0; ki<K; ki++)
			    {
				Di[ki] += lambda*m_V[kj*K + ki];
			    }
		    }




		if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
		if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

		// compute min
		delta = Di[0];
		m_answer[n] = 0;
		for (ki=1; ki<K; ki++)
		    {
			if (delta > Di[ki])
			    {
				delta = Di[ki];
				m_answer[n] = ki;
			    }
		    }
	    }

    delete [] Di;
}

void BPS::computeSolution_GRID_BINARY() {

    int x, y, n;
    REAL* M_ptr;
    REAL Di;

    n = 0;
    M_ptr = m_messages;

    for (y=0; y<m_height; y++)
	for (x=0; x<m_width; x++, M_ptr+=2, n++)
	    {
		Di = m_DBinary[n];
		if (x > 0) Di += (m_answer[n-1] == 0)       ? m_horzWeightsBinary[n-1]       : -m_horzWeightsBinary[n-1];
		if (y > 0) Di += (m_answer[n-m_width] == 0) ? m_vertWeightsBinary[n-m_width] : -m_vertWeightsBinary[n-m_width];

		if (x < m_width-1)  Di += M_ptr[0]; // message (x+1,y)->(x,y)
		if (y < m_height-1) Di += M_ptr[1]; // message (x,y+1)->(x,y)

		// compute min
		m_answer[n] = (Di >= 0) ? 0 : 1;
	    }

}

void BPS::computeSolution_GRID_GENERAL() {

    int x, y, n, K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    REAL delta;
    int ki, kj;

    Di = new REAL[K];

    n = 0;
    D_ptr = m_D;
    M_ptr = m_messages;

    for (y=0; y<m_height; y++)
	for (x=0; x<m_width; x++, D_ptr+=K, M_ptr+=2*K, n++)
	    {
		CopyVector(Di, D_ptr, K);
		if (m_V)
		    {
			CostVal* ptr = m_V + 2*(x+y*m_width-1)*K*K;
			if (x > 0)
			    {
				kj = m_answer[n-1];
				for (ki=0; ki<K; ki++)
				    {
					Di[ki] += ptr[kj + ki*K];
				    }
			    }
			ptr -= (2*m_width-3)*K*K;
			if (y > 0)
			    {
				kj = m_answer[n-m_width];
				for (ki=0; ki<K; ki++)
				    {
					Di[ki] += ptr[kj + ki*K];
				    }
			    }
		    }
		else
		    {
			if (x > 0)
			    {
				kj = m_answer[n-1];
				for (ki=0; ki<K; ki++)
				    {
					Di[ki] += m_smoothFn(n, n-1, ki, kj);
				    }
			    }
			if (y > 0)
			    {
				kj = m_answer[n-m_width];
				for (ki=0; ki<K; ki++)
				    {
					Di[ki] += m_smoothFn(n, n-m_width, ki, kj);
				    }
			    }
		    }

		if (x < m_width-1) AddVector(Di, M_ptr, K); // message (x+1,y)->(x,y)
		if (y < m_height-1) AddVector(Di, M_ptr+K, K); // message (x,y+1)->(x,y)

		// compute min
		delta = Di[0];
		m_answer[n] = 0;
		for (ki=1; ki<K; ki++)
		    {
			if (delta > Di[ki])
			    {
				delta = Di[ki];
				m_answer[n] = ki;
			    }
		    }
	    }

    delete [] Di;


}

void BPS::computeSolution_NONGRID_FIXED_MATRIX_GENERAL() {

    int K = m_nLabels;
    CostVal* D_ptr;
    REAL* M_ptr;
    REAL* Di;
    REAL delta;
    int ki, kj;

    CostVal lambda = 1;

    Di = new REAL[K];

    D_ptr = m_D;

    for (int m=0; m<m_nPixels; m++, D_ptr+=K) {
    	CopyVector(Di, D_ptr, K);

    	// for the neighbors that are already assigned labels use those labels
    	// this check is performed by seeing if label = K
    	// which means that node is not yet labeled
    	// for neighbors not assigned labels, use the messages instead

    	M_ptr = m_messagesNG[m];

    	for (int p=0; p<m_nNeighbors[m]; p++, M_ptr+=K) {

    		int q=findNeighborWithIdx(m, p);
    		if (q==-1) {
    			printf("Error\n");
    			exit(1);
    		}

    		if (m_answer[q] == K) { //neighbor not labeled
    			AddVector(Di, M_ptr, K);
    		} else { // neighbor is labeled
    			if ( m_smoothType != FUNCTION  ) { //assume ARRAY
					kj = m_answer[q];
					for (ki=0; ki<K; ki++)
					{
						Di[ki] += lambda*m_V[kj*K + ki];
					}
    			} else { // function
					kj = m_answer[q];
					for (ki=0; ki<K; ki++)
					{
						Di[ki] += lambda*m_smoothFn(m, q, ki, kj);
					}

    			}
    		}
    	}

    	//then compute min

		delta = Di[0];
		m_answer[m] = 0;
		for (ki=1; ki<K; ki++)
		{
			if (delta > Di[ki])
			{
				delta = Di[ki];
				m_answer[m] = ki;
			}
	    }
//		printf("Answer pix %d = %d\n",m, m_answer[m]);
    }

    delete [] Di;

}


void BPS::fillAnswersK() {
	for (int i=0; i<m_nPixels; i++) {
		m_answer[i] = m_nLabels;
	}
}

