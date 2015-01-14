#include "LinkedBlockList.h"
#include<iostream>

/*********************************************************************/


void LinkedBlockList::addFront(ListType item) {

    if ( m_head_block_size == GCLL_BLOCK_SIZE )
    {
        LLBlock *tmp      = (LLBlock *) new LLBlock;
        tmp -> m_next     = m_head;
        m_head            = tmp;
        m_head_block_size = 0;
    }

    m_head ->m_item[m_head_block_size] = item;
    m_head_block_size++;
}

/*********************************************************************/

void LinkedBlockList::addFront(ListType item, int mflag) {

	m_neighborFlag = mflag; // indicating that data type of item is of class Neighbor;
	// temporary bug fix //chet

    if ( m_head_block_size == GCLL_BLOCK_SIZE )
    {
        LLBlock *tmp      = (LLBlock *) new LLBlock;
        tmp -> m_next     = m_head;
        m_head            = tmp;
        m_head_block_size = 0;
    }

    m_head ->m_item[m_head_block_size] = item;
    m_head_block_size++;
}


/*********************************************************************/

ListType LinkedBlockList::next()
{
    ListType toReturn = m_cursor -> m_item[m_cursor_ind];

    m_cursor_ind++;

    if ( m_cursor == m_head && m_cursor_ind >= m_head_block_size )
    {
        m_cursor     = m_cursor ->m_next;
        m_cursor_ind = 0;
    }
    else if ( m_cursor_ind == GCLL_BLOCK_SIZE )
    {
        m_cursor = m_cursor ->m_next;
        m_cursor_ind = 0;
    }
    return(toReturn);
}

/*********************************************************************/

bool LinkedBlockList::hasNext()
{
    if ( m_cursor != 0 ) return (true);
    else return(false);
}


/*********************************************************************/

//chet
// these are the definitions in successive parent files for defining the class Neighbor
/*
typedef float CostVal;
typedef MRF::CostVal captype;
typedef Graph::captype EnergyTermType;
typedef int PixelType;
typedef struct NeighborStruct{
	PixelType  to_node;
	EnergyTermType weight;
} Neighbor;
*/

LinkedBlockList::~LinkedBlockList()
{

	typedef int PixelType;
	typedef float EnergyTermType;
	typedef struct NeighborStruct{
		PixelType  to_node;
		EnergyTermType weight;
	} Neighbor;

    LLBlock *tmp;

    int cnt = 0;

    while ( m_head != 0 )
    {
//chet
    	if (cnt==0) {
    		for (int i=0; i<m_head_block_size; i++) {
				if (m_neighborFlag==1)
					delete (Neighbor*) m_head ->m_item[i];
				else {
					delete m_head ->m_item[i];
					std::cout << "C++ does not define deleting void ptrs, so could be source of segmentation fault. Not including this line will cause memory leaks \n";
				}
			}
    	} else {
    		for (int i=0; i<GCLL_BLOCK_SIZE; i++) {
				if (m_neighborFlag==1)
					delete (Neighbor*) m_head ->m_item[i];
				else {
					delete m_head ->m_item[i];
					std::cout << "C++ does not define deleting void ptrs, so could be source of segmentation fault. Not including this line will cause memory leaks \n";
				}
			}
		}
    	cnt++;
        tmp = m_head;
        m_head = m_head->m_next;
        delete tmp;
    }
};

/*********************************************************************/

