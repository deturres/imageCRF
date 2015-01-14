#include "imageLib.h"
#include "mrf.h"
#include "LabelImageWriter.h"


LabelImageWriter::LabelImageWriter(char* stem, int width, int height, float scale) : 
  m_iter(0) 
  {
    m_stem = stem;
    m_width = width;
    m_height = height;
    m_scale = scale;
  }

void LabelImageWriter::write(MRF::Label* labels) {

  printf("writing image\n");

  CByteImage im;
  
  CShape sh(m_width, m_height, 1);
  im.ReAllocate(sh);
  
  for (int y = 0; y < m_height; y++) 
    {
      uchar *row = &im.Pixel(0, y, 0);
      for (int x = 0; x < m_width; x++, labels++) 
        {
          row[x] = *labels;
        }
    }
  
  if (m_scale!=1)
    ScaleAndOffset(im,im,m_scale,0);
  
  char name[1000];
  sprintf(name,"%s-%d.png",m_stem,m_iter);
  WriteImage(im,name);
  m_iter++;
}
