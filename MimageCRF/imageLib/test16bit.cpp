//
//  test16bit.cpp
//  
//
//  Created by Chetan on 12/5/12.
//
//

#include "test16bit.h"

#include <iostream>
using namespace std;

#include "imageLib.h"
#include "Image.h"

int main() {
  
  CImage2 img2;
  char fname[100] = "im16bit.png";
  read_png_file(img2, fname);
  
  printf("img2 200, 235, 0 = <%d>\n", (int)img2.Pixel(200,235,0));
	printf("img2 200, 235, 1 = <%d>\n", (int)img2.Pixel(200,235,1));
	printf("img2 200, 235, 2 = <%d>\n", (int)img2.Pixel(200,235,2));
	printf("img2 201, 235, 0 = <%d>\n", (int)img2.Pixel(201,235,0));
	printf("img2 201, 235, 1 = <%d>\n", (int)img2.Pixel(201,235,1));
	printf("img2 201, 235, 2 = <%d>\n", (int)img2.Pixel(201,235,2));
	printf("img2 200, 236, 0 = <%d>\n", (int)img2.Pixel(200,236,0));
	printf("img2 200, 236, 1 = <%d>\n", (int)img2.Pixel(200,236,1));
	printf("img2 200, 236, 2 = <%d>\n", (int)img2.Pixel(200,236,2));

  CImage2 img1;
  char fname2[100] = "imgray16bit.png";
  read_png_file(img1, fname2);
  
  printf("img1 200, 235, 0 = <%d>\n", (int)img1.Pixel(200,235,0));
	printf("img1 201, 235, 0 = <%d>\n", (int)img1.Pixel(201,235,0));
	printf("img1 200, 236, 0 = <%d>\n", (int)img1.Pixel(200,236,0));
  
  printf("img1 639, 479, 0 = <%d>\n", (int)img1.Pixel(639,479,0));

  return 0;
}
