#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mrf.h"

double DataCost::getDataCostRaw(int index)
{
	if (m_type==MRF::ARRAY)
	{
		return (double) m_costArray[index];
	} else
	{
		exit(-1);
	}
}

int DataCost::updateDataCostRaw(int index, double newvalue)
{
	if (m_type==MRF::ARRAY)
	{
		m_costArray[index] = newvalue ;
	} else
	{
		exit(-1);
	}

}


double SmoothnessCost::getSmoothnessCostRaw(int p, int q, int lp, int lq)
{
	if (m_type == MRF::FUNCTION)
	{
		return (double) (*m_costFn) (p, q, lp, lq);
	} else
	{
		exit(-1);
	}

}


void MRF::initialize()
{

    m_initialized = 1;

    if ( m_dataType == ARRAY )
        setData(m_e->m_dataCost->m_costArray);
    else  setData(m_e->m_dataCost->m_costFn);

    if ( m_smoothType == FUNCTION )
        setSmoothness(m_e->m_smoothCost->m_costFn);
    else 
    {
        if ( m_smoothType == ARRAY )
        {
            checkArray(m_e->m_smoothCost->m_V);
            setSmoothness(m_e->m_smoothCost->m_V);
        }
        else
        {
            int smoothExp = m_e->m_smoothCost->m_smoothExp;
            if ( smoothExp != 1 && smoothExp != 2) 
            { 
                fprintf(stderr, "Wrong exponent in setSmoothness()!\n"); 
                exit(1); 
            }
            setSmoothness(smoothExp,m_e->m_smoothCost->m_smoothMax,m_e->m_smoothCost->m_lambda);
        }

        if (m_e->m_smoothCost->m_varWeights )
        {
            if (!m_grid_graph) 
            { 
                fprintf(stderr, "Edge multiplier cannot be used with non-grid graphs\n"); 
                exit(1); 
            }
            setCues(m_e->m_smoothCost->m_hWeights,m_e->m_smoothCost->m_vWeights);
        }

    }

	initializeAlg();
}


void MRF::commonInitialization(EnergyFunction *e)
{
    m_dataType    = e->m_dataCost->m_type;
    m_smoothType  = e->m_smoothCost->m_type;
    m_varWeights  = e->m_smoothCost->m_varWeights;
    m_initialized = 0;
    m_e = e;
    m_allocateArrayForSmoothnessCostFn = true;
}


MRF::MRF(int width, int height, int nLabels, EnergyFunction *e)
{
    m_width       = width;
    m_height      = height;
    m_nLabels     = nLabels;
    m_nPixels     = width*height;
    m_grid_graph  = 1;
    commonInitialization(e);
}

MRF::MRF(int nPixels, int nLabels, EnergyFunction *e)
{
    m_nLabels     = nLabels;
    m_nPixels     = nPixels;
    m_grid_graph  = 0;
    commonInitialization(e);
}

char MRF::checkEnergy()
{ 
    if ( !m_initialized) {fprintf(stderr, "Call initialize() first,exiting!");exit(1);}
    return(2);
}


MRF::EnergyVal MRF::totalEnergy()
{
    if (!isValid()) { fprintf(stderr, "totalEnergy() cannot be called for invalid energy!\n"); exit(1); }
    if (!m_initialized) { fprintf(stderr, "Call initialize() first!\n"); exit(1); }

    return dataEnergy()+smoothnessEnergy();
}


void MRF::optimize(int nIterations, float& time)
{
    if (!isValid()) { fprintf(stderr, "optimize() cannot be called for invalid energy!\n"); exit(1); }
    if (!m_initialized ) { fprintf(stderr, "run initialize() first!\n"); exit(1); }

    // start timer
    start = clock();
    timeoffset = 0;

    optimizeAlg(nIterations);

    // stop timer
    clock_t finish = clock();
    time = (float) (((double)(finish - start + timeoffset)) / CLOCKS_PER_SEC);
}


void MRF::checkArray(CostVal *V)
{

    int i, j;
    
    for (i=0; i< m_nLabels; i++)
        for (j=i; j<m_nLabels; j++)
        if (V[i*m_nLabels + j] != V[j*m_nLabels + i])
        { 
            fprintf(stderr, "Error in setSmoothness(): V is not symmetric!\n"); 
            exit(1); 
        }
        

}
