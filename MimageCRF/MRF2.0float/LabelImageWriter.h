#ifndef __LABELIMAGEWRITER_H__
#define __LABELIMAGEWRITER_H__

class LabelImageWriter {

public:

  LabelImageWriter(char* stem, int width, int height, float scale=1);

  void write(MRF::Label* labels);
  
  void setStem(char* stem)  { m_stem = stem; };

 private:
  char* m_stem;
  int m_iter, m_width, m_height;
  float m_scale;
};

#endif
