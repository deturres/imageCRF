/*
 * IO.h
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#ifndef IO_H_
#define IO_H_

#include <string>
#include <vector>
#include "imageLib/Error.h"
#include <stdio.h>
#include "helperLib.h"
#include "Parameters.h"

class IO
{

	int copyDirList(std::vector<string> dirList);


protected:

	char debugname[500];
	// char logname[500];

	int numTrainingPats;
	int numTestPats;

	int timeseq;

	std::vector<string> dirListNFP; //input patient folder names (no path)
  std::vector<string> indirListFP; //input patient folder names (fullpath)
  std::vector<string> gtdirListFP; //groundtruth patient folder names (fullpath)
  std::vector<string> interdirListFP; //prm.interactive map patient folder names (fullpath)
  std::vector<string> outdirListFP; //output patient folder names (fullpath) - need to create these manually for now.

  std::vector <std::vector <CImage2> > indirImage; // contains all images of all patients and all slices
  std::vector <std::vector <CByteImage> > indirImaget; // contains all images of all videos and all frames
  std::vector <std::vector <CByteImage> > gtdirImage; // contains gtruth of all patients and all slices or videos and frames
  std::vector <std::vector <CByteImage> > interdirImage; // contains prm.interactive map of all patients and all slices or videos and all frames

  std::vector <std::vector <CByteImage> > dirDisp; //images of output for all patients and all slices
  std::vector <std::vector <CByteImage> > dirDispGen; //image with prm.featurepv.location enabled output

  std::vector <std::vector <CByteImage> > errormap; //images of output for all patients and all slices
  std::vector <std::vector <CByteImage> > errormapGen; //image with prm.featurepv.location enabled output

public:

	FILE *debugfile;
	FILE *paramoutfile;
	FILE *resultoutfile;

	IO(std::string outstem, int tseq);
	~IO();
	int readIOfolderStructures(ioparams &iopv, int interactive, int crfhidden);
	int readIOfiles(int interactive, int crfhidden);
	int deAllocateIOfiles();
	int allocateOutputSpace(int generative);

	int getNumPats()
	{
		return numTrainingPats+numTestPats;
	}

	std::vector <std::vector <CImage2> > * getInputImageDirRef();
	std::vector <std::vector <CByteImage> > * getInputImagetDirRef();
	std::vector <std::vector <CByteImage> > * getGtImageDirRef();
	std::vector <std::string> getVolumeNames();
	std::vector <string> getInputImagetDirPath();

	int scaleGtruth(int outscale8);
	int scaleGtruthWithUnknown(int outscale8, int nD);

};




#endif /* IO_H_ */
