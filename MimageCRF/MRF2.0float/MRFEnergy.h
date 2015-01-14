/**
 * Class extending MRF but containing all the standard useful methods for
 * dealing with the energy function.
 *
 * Jerod Weinman
 * jerod@acm.org
 */


#ifndef __MRFENERGY_H__
#define __MRFENERGY_H__

#include "mrf.h"
#include "LabelImageWriter.h"
#include "LinkedBlockList.h"
#include <ostream>

#define FloatType float

class MRFEnergy : public MRF {

 public:

  // Function handle for writing images
  typedef void (*ImageWriteFn)(Label* labels);

  MRFEnergy(int width, int height, int nLabels, EnergyFunction *eng);

  MRFEnergy(int nPixels, int nLabels,EnergyFunction *eng);

  ~MRFEnergy();

  // Returns the total energy either for smoothness or data
  EnergyVal smoothnessEnergy();
  EnergyVal dataEnergy();

  // Data cost from the energy function for a particular assignment
  CostVal getDataCost(const int pixel, const Label label);

  // Smoothness cost from the energy function for a particular assignment
  CostVal getSmoothnessCost(int pixel1, int pixel2,
                            Label label1, Label label2);

  // Copies the energy terms into a destination vector
  void copyDataCost(int pixel, CostVal* dst);

  // Adds the smoothness costs to a destination vector or matrix
  void addSmoothnessCost(const int fixedPixel, const int varPixel,
                         const Label fixedLabel,
                         CostVal* dst);

  void addSmoothnessCost(const int pixel1, const int pixel2, CostVal* dst);

  // Copies the energy terms into a destination vector in NON-SPARSE positions
  void copySparseDataCost(int pixel, CostVal* dst, bool* sp);

  // Adds the smoothness costs to a destination vector in NON-SPARSE positions
  void addSparseSmoothnessCost(const int fixedPixel, const int varPixel,
                               const Label fixedLabel,
                               CostVal* dst, bool *sp);


  // Indicates that spatially varying weights on the smoothness term
  // (edge potential) are being used
  bool varWeights();

  // Indicates whether a function, lookup table, etc. is used for smoothness
  InputType getSmoothType();

  // Fetch values of spatially varying weights on the
  // smoothness matrix (edge potentials)
  CostVal getHorizWeight(const int r, const int c);
  CostVal getVertWeight(const int r, const int c);

  // Returns pointer to space to write a "solution"
  Label* getAnswerPtr();


  Label getLabel(int pixel);
  void setLabel(int pixel,Label label);
  void clearAnswer();

  inline int getNumPixels() { return m_nPixels; }
  inline int getNumLabels() { return m_nLabels; }

  // Get the current approximate marginal for a pixel *in logspace*
  virtual void getBelief(const int pixel, FloatType* dst)=0;

  virtual FloatType getBelief(const int pixel, const Label label)=0;

  // Get the current approximate marginal for two neighboring pixels *in logspace*
  virtual void getEdgeBelief(const int pixel1, const int pixel2, FloatType* dst)=0;

  // Get the log probability of a configuration
  virtual FloatType getLogProb(Label* labels)=0;

  // Linear index from 2-D subscripts (row-major order)
  inline int pixelIndex(int row, int col) { return row*m_width + col; };

  // Linear pairwise label index from 2-D subscripts
  inline int labelIndex(Label label1, Label label2)
    { return label1*m_nLabels + label2; }

#ifdef ENERGY_FRIEND_COSTS // This condition shouldn't be met probably
  void setData(DataCost* dcost);
  void setSmoothness(SmoothnessCost*);
#endif

  // Option in case derived classes want to write intermediate output images
  void setLabelImageWriter(LabelImageWriter *imWriter);

  // Publicize these so we can update the energy function
  // (in a learning framework)
  void setData(DataCostFn dcost);
  void setData(CostVal* data);
  void setSmoothness(SmoothCostGeneralFn cost);
  void setSmoothness(CostVal* V);
  void setSmoothness(int smoothExp,CostVal smoothMax, CostVal lambda);
  void setCues(CostVal* hCue, CostVal* vCue);

  // Set an output stream for logging purposes
  void setLog(std::ostream *stream);

  virtual void initializeAlg()=0;

//chet
  virtual void setLog2(std::ostream *stream) { }
  void setNeighbors(int pix1, int pix2, CostVal weight);
  int traverseNeighbors(int pixel1);


 protected:

  CostVal *m_V; // Cost array for smoothness terms
  CostVal *m_D; // Cost array for data terms
  CostVal *m_horizWeights;
  CostVal *m_vertWeights;
  bool m_varWeights;
  DataCostFn m_dataFn;
  SmoothCostGeneralFn m_smoothFn;
  LabelImageWriter *m_imWriter; // For writing intermediate results

  // Prediction resulting from running the algorithm
  Label* m_labels;

  // Indicates the smoothness cost array was allocated (and thus needs deleting)
  bool m_freeV;

  // Log stream to write output to
  std::ostream *m_log;

//chet
	typedef struct NeighborStruct{
	 int  to_node;
	 CostVal weight;
	} Neighbor;
	LinkedBlockList *m_neighbors;

	int *m_nNeighbors;
	int m_maxnNB;

    int findIdxOfNeighbor(int mainpix, int pixel1);
    int findNeighborWithIdx(int mainpix, int idx);
	void fillAnswersK();



 private:

  // Common initialization routine
  void initMRFEnergy();


};


#endif
