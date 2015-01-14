#include <iostream>
using namespace std;

#include "imageLib.h"
#include "Image.h"

//extern "C"{
//#include "png.h"
//}

//#include "png.h"
//#include "ImageIO.h"
//#include "Convert.h"
//#include "Error.h"



int main()
{

	CImage2 img;
	CImage2 img2;
	CByteImage img3;

	read_png_file(img, "16_1.png");
	read_png_file(img2, "16_3.png");

	ReadImage(img3, "seg.png");


	printf("img_1 200, 235, 0 = <%d>\n", (int)img.Pixel(200,235,0));
	printf("img_1 201, 235, 0 = <%d>\n", (int)img.Pixel(201,235,0));
	printf("img_1 200, 236, 0 = <%d>\n", (int)img.Pixel(200,236,0));	
	printf("img_1 235, 200, 0 = <%d>\n", (int)img.Pixel(235,200,0));
	
	printf("img_3 200, 235, 0 = <%d>\n", (int)img2.Pixel(200,235,0));	
	printf("img_3 200, 235, 1 = <%d>\n", (int)img2.Pixel(200,235,1));	
	printf("img_3 200, 235, 2 = <%d>\n", (int)img2.Pixel(200,235,2));
	printf("img_3 201, 235, 0 = <%d>\n", (int)img2.Pixel(201,235,0));	
	printf("img_3 201, 235, 1 = <%d>\n", (int)img2.Pixel(201,235,1));	
	printf("img_3 201, 235, 2 = <%d>\n", (int)img2.Pixel(201,235,2));
	printf("img_3 200, 236, 0 = <%d>\n", (int)img2.Pixel(200,236,0));	
	printf("img_3 200, 236, 1 = <%d>\n", (int)img2.Pixel(200,236,1));	
	printf("img_3 200, 236, 2 = <%d>\n", (int)img2.Pixel(200,236,2));
	
	printf("img_3 235, 200, 0 = <%d>\n", (int)img2.Pixel(235,200,0));	
	printf("img_3 235, 200, 1 = <%d>\n", (int)img2.Pixel(235,200,1));	
	printf("img_3 235, 200, 2 = <%d>\n", (int)img2.Pixel(235,200,2));
	
	printf("seg 200, 235, 0 = <%d>\n", (int)img3.Pixel(200,235,0));	
	printf("seg 235, 200, 0 = <%d>\n", (int)img3.Pixel(235,200,0));



	return 0;
}
