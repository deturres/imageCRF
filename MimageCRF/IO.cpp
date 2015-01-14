/*
 * IO.cpp
 *
 *  Created on: Mar 15, 2011
 *      Author: bhole
 */

#include "IO.h"


IO::IO(std::string outstem, int tseq)
{

  // debug file
	sprintf(debugname, "%s.txt", outstem.c_str());
	debugfile = fopen(debugname, "w");

	if (debugfile == NULL)
	  throw CError("Cannot write to %s\n", debugname);

	char filename[500];
	// output param file
	sprintf(filename, "%sparam.txt", outstem.c_str());
	paramoutfile = fopen(filename, "w");

	if (paramoutfile == NULL)
	  throw CError("Cannot write to %s\n", filename);


	// output result file
	sprintf(filename, "%sresult.txt", outstem.c_str());
	resultoutfile = fopen(filename, "w");

	if (resultoutfile == NULL)
	  throw CError("Cannot write to %s\n", filename);



	numTrainingPats = 0;
	numTestPats = 0;

	timeseq = tseq;

}

IO::~IO()
{
	fclose(debugfile);
	fclose(paramoutfile);
	fclose(resultoutfile);
}


int IO::copyDirList(std::vector<string> dirList)
{
	dirListNFP = dirList;
}

std::vector <std::string> IO::getVolumeNames()
{
	return dirListNFP;
}


int IO::readIOfolderStructures(ioparams &iopv, int interactive, int crfhidden)
{

    std::vector<string> dirList; // temporary variable: stores the subdirectory names (patient folder names)

    // input
    dirList.erase(dirList.begin(), dirList.end());
    readDirectories(iopv.indirname, dirList, indirListFP, iopv.testDir);

    // ground-truth
    dirList.erase(dirList.begin(), dirList.end());
    readDirectories(iopv.gtdirname, dirList, gtdirListFP, iopv.testDir);

    copyDirList(dirList);

    if (interactive == 1 || crfhidden == 1)
    {
    	dirList.erase(dirList.begin(), dirList.end());
    	readDirectories(iopv.interdirname, dirList, interdirListFP, iopv.testDir);
    }

    // output
    dirList.erase(dirList.begin(), dirList.end());
    readDirectories(iopv.outdirname, dirList, outdirListFP, iopv.testDir);

    iopv.testDirIndexV.resize(dirList.size());

    // find out which are test and which are training
    if (iopv.testDir == 1)
    {
      int jjj=0;
      for(std::vector<string>::iterator b = dirList.begin(); b != dirList.end(); ++b)
      {
        string foldn = *b;
        if (foldn.find("test") == 0)
        {
          iopv.testDirIndexV[jjj] = 1;
          numTestPats++ ;
        }
        else {
          iopv.testDirIndexV[jjj] = 0;
          numTrainingPats++ ;
        }

        jjj++;
      }
    } else {

      numTrainingPats = dirList.size();

      int jjj=0;
      for(std::vector<string>::iterator b = dirList.begin(); b != dirList.end(); ++b)
      {
        iopv.testDirIndexV[jjj] = 0;
        jjj++;
      }
    }


    // check for consistency

    // the number of sub-directories (patient folders) must match
    if (indirListFP.size() != gtdirListFP.size())
      throw CError("Not all number of directories are the same (Check input and gtruth directories) .\n");

    if (indirListFP.size() != outdirListFP.size())
      throw CError("Not all number of directories are the same (Check input and output directories) .\n");

    // the number of sub-directories (patient folders) must match
    if (interactive == 1 || crfhidden == 1)
        if( (indirListFP.size() != interdirListFP.size()))
          throw CError("Not all number of directories are the same (Check interactive directories).\n");

    if (iopv.testDir==1 && iopv.testDirIndexV.size()==0)
      throw CError("test directory missing or parameter prm.iopv.testDir falsely enabled.\n");


    return 0;

}


int IO::deAllocateIOfiles()
{
    if (timeseq == 0)
    {
    	int numD = indirImage.size();

    	for (int i=0; i<numD; i++)
    	{
    		int numD2 = indirImage[i].size();
    		for (int j=0; j<numD2; j++)
    		{
    			// std::cout << "Deallocating " <<i <<"  " << j <<"\n";
    			indirImage[i][j].DeAllocate();
    		}
    	}
    }

}



int IO::readIOfiles(int interactive, int crfhidden)
{

  std::vector<CImage2> set1Image;
  std::vector<CByteImage> set1Imaget;
  std::vector<CByteImage> setTrueDisp;
  std::vector<CByteImage> setInterDisp;

  indirImage.erase(indirImage.begin(), indirImage.end());
  indirImaget.erase(indirImaget.begin(), indirImaget.end());

  for(std::vector<string>::iterator b = indirListFP.begin(); b < indirListFP.end(); ++b)
  {
    std::cout << *b <<"\n";
    set1Image.erase(set1Image.begin(), set1Image.end());
    set1Imaget.erase(set1Imaget.begin(), set1Imaget.end());
    if (timeseq == 0)
    {
      readImages(*b, set1Image);
      indirImage.push_back(set1Image);
    } else if (timeseq == 1) {
      readImages(*b, set1Imaget);
      indirImaget.push_back(set1Imaget);
    }
  }

  for(std::vector<string>::iterator b = gtdirListFP.begin(); b < gtdirListFP.end(); ++b)
  {
    cout << *b <<"\n";
    setTrueDisp.erase(setTrueDisp.begin(), setTrueDisp.end());
    readImages(*b, setTrueDisp);
    gtdirImage.push_back(setTrueDisp);
  }

  if (interactive == 1 || crfhidden == 1)
  {
    for(std::vector<string>::iterator b = interdirListFP.begin(); b < interdirListFP.end(); ++b)
    {
      cout << *b <<"\n";
      setInterDisp.erase(setInterDisp.begin(), setInterDisp.end());
      readImages(*b, setInterDisp);
      interdirImage.push_back(setInterDisp);
    }
  } else
  {
    for(std::vector<string>::iterator b = gtdirListFP.begin(); b < gtdirListFP.end(); ++b)
      interdirImage.push_back((std::vector<CByteImage>) 0);
  }

  // Check folders correspond and also need to check consistency of number of images in all
  // folders of corresponding input/gtruth and/or interactive should be same

  if (timeseq == 0)
  {
    assert(gtdirImage.size() == indirImage.size());
    for (unsigned int i = 0; i < gtdirImage.size(); i++)
    {
      assert(gtdirImage[i].size() == indirImage[i].size());
    }
  }
  if (timeseq == 1)
  {
    assert(gtdirImage.size() == indirImaget.size());
    for (unsigned int i = 0; i < gtdirImage.size(); i++)
    {
      assert(gtdirImage[i].size() == indirImaget[i].size());
    }
  }

  assert(gtdirImage.size() == interdirImage.size()); // forced even when interactive = 0;
  if (interactive == 1 || crfhidden == 1)
  {
    for (unsigned int i = 0; i < gtdirImage.size(); i++)
    {
      assert(gtdirImage[i].size() == interdirImage[i].size());
    }
  }



	return 0;
}

int IO::scaleGtruth(int outscale8)
{
    for(unsigned int j = 0; j < gtdirImage.size(); j++) {
      for(unsigned int i = 0; i < gtdirImage[j].size(); ++i) {
        ScaleAndOffset(gtdirImage[j][i], gtdirImage[j][i], (double)1.0/outscale8, 0); // scale to integer disps
      }
    }
}

// This function will handle the unknown case.
int IO::scaleGtruthWithUnknown(int outscale8, int nD)
{
    for(unsigned int j = 0; j < gtdirImage.size(); j++) {
      for(unsigned int i = 0; i < gtdirImage[j].size(); ++i) {
        ScaleAndOffsetWithUnknown(gtdirImage[j][i], gtdirImage[j][i], (double)1.0/outscale8, 0, nD); // scale to integer disps
      }
    }
}


int IO::allocateOutputSpace(int generative)
{

  dirDisp.resize(gtdirImage.size());
  errormap.resize(gtdirImage.size());

  if (generative==1)
  {
    dirDispGen.resize(gtdirImage.size());
    errormapGen.resize(gtdirImage.size());
  }

  for (unsigned int j = 0; j < gtdirImage.size(); ++j)
  {

    int depth = gtdirImage[j].size();
    CShape sh = gtdirImage[j][0].Shape();
    int width = sh.width, height = sh.height, nB = sh.nBands;
    sh.nBands = 1;

    dirDisp[j].resize(gtdirImage[j].size());
    if (generative==1)
    {
      dirDispGen[j].resize(gtdirImage[j].size());
    }

    for(unsigned int i = 0; i < depth; ++i)
    {
      dirDisp[j][i].ReAllocate(sh);

      // set all values to 0
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          uchar *row = &dirDisp[j][i].Pixel(x, y, 0);
          *row = 0;
        }
      }

      if (generative==1)
      {
        dirDispGen[j][i].ReAllocate(sh);

        // set all values to 0
        for (int y = 0; y < height; y++)
        {
          for (int x = 0; x < width; x++)
          {
            uchar *row = &dirDispGen[j][i].Pixel(x, y, 0);
            *row = 0;
          }
        }
      }
    }
  }
}




std::vector <std::vector <CImage2> > * IO::getInputImageDirRef()
{
	return &indirImage;
}

std::vector <std::vector <CByteImage> > * IO::getInputImagetDirRef()
{
	return &indirImaget;
}

std::vector <string> IO::getInputImagetDirPath()
{
  return indirListFP;
}


std::vector <std::vector <CByteImage> > * IO::getGtImageDirRef()
{
	return &gtdirImage;
}
