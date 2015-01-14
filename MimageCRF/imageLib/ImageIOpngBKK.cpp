///////////////////////////////////////////////////////////////////////////
//
// NAME
//  ImageIOpng.cpp -- reads and writes png images
//
// DESCRIPTION
//  This file is based pmon lpng\contrib\visupng\PngFile.c by Willem van Schaik
//  It requires the libpng and zlib libraries.
//
// SEE ALSO
//  ImageIO.cpp
//
// Copyright � Richard Szeliski and Daniel Scharstein, 2005.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

extern "C"{
#include "png.h"
}

#include <stdarg.h>
#include "Image.h"
#include "Error.h"
#include <vector>

// global variables pointing to the png stuctures
static png_structp png_ptr = NULL;
static png_infop info_ptr = NULL;


// chetan
// both these functions are the same since it just means interchanging the 2 bytes.
unsigned short endianl2b(unsigned short i)
{
    return ((i>>8)&0xff)+((i << 8)&0xff00);
}

unsigned short endianb2l(unsigned short i)
{
    return ((i>>8)&0xff)+((i << 8)&0xff00);
}

void abort_(const char * s, ...)
{
	va_list args;
	va_start(args, s);
	vfprintf(stderr, s, args);
	fprintf(stderr, "\n");
	va_end(args);
	abort();
}















static void pngfile_error(png_structp /*png_ptr*/, png_const_charp msg)
{
	throw CError(msg);
}

#define DEBUG_ImageIOpng 0

// TODO: the following function should go somewhere else, perhaps in ImageIO.cpp
// Make sure the image has the smallest number of bands before writing.
// That is, if it's 4 bands with full alpha, reduce to 3 bands.  
// If it's 3 bands with constant colors, make it 1-band.
CByteImage removeRedundantBands(CByteImage img)
{
    CShape sh = img.Shape();
	int w = sh.width, h = sh.height, nB = sh.nBands;
	int x, y;
	if (nB < 3)
		return img;

	// check if full alpha if alpha channel present
	bool fullAlpha = true;
	if (nB == 4) {
		for (y = 0; y < h && fullAlpha; y++) {
			uchar *pix = &img.Pixel(0, y, 0);
			for (x = 0; x < w; x++) {
				if (pix[3] != 255) {
					fullAlpha = false;
					break;
				}
				pix += nB;
			}
		}
	}
	if (!fullAlpha)
		return img;

	// check for equal colors
	bool equalColors = true;
	for (y = 0; y < h && equalColors; y++) {
		uchar *pix = &img.Pixel(0, y, 0);
		for (x = 0; x < w; x++) {
			if (pix[0] != pix[1] ||
				pix[0] != pix[2] ||
				pix[1] != pix[2]) {
					equalColors = false;
					break;
				}
				pix += nB;
		}
	}
	// at this point, if nB == 4 we can reduce to at least 3 bands,
	// and if equalColors we can reduce to 1 band.
	if (! equalColors && nB < 4)
		return img;

	int newNB = equalColors ? 1 : 3;

	if (DEBUG_ImageIOpng)
		fprintf(stderr, "reducing from %d to %d bands\n", nB, newNB);

	CShape sh2(w, h, newNB);
	CByteImage img2(sh2);
	
	for (y = 0; y < h; y++) {
		uchar *pix = &img.Pixel(0, y, 0);
		uchar *pix2 = &img2.Pixel(0, y, 0);
		for (x = 0; x < w; x++) {
			for (int b = 0; b < newNB; b++) {
				pix2[b] = pix[b];
			}
			pix += nB;
			pix2 += newNB;
		}
	}

	return img2;
}


//chetan

void read_png_file(CImage2& img, char* file_name)
{
	int x,y;

	int width, height;
	png_byte color_type;
	png_byte bit_depth;

	//png_structp png_ptr;
	//png_infop info_ptr;
	int number_of_passes;
	png_bytep * row_pointers;

	char header[8];	// 8 is the maximum size that can be checked

	/* open file and test for it being a png */
	FILE *fp = fopen(file_name, "rb");
	if (!fp)
		abort_("[read_png_file] File %s could not be opened for reading", file_name);
	fread(header, 1, 8, fp);
	if (png_sig_cmp((png_byte *) header, 0, 8))
		abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


	/* initialize stuff */
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	
	if (!png_ptr)
		abort_("[read_png_file] png_create_read_struct failed");

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		abort_("[read_png_file] png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[read_png_file] Error during init_io");

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	width = info_ptr->width;
	height = info_ptr->height;
	color_type = info_ptr->color_type;
	bit_depth = info_ptr->bit_depth;

	//printf("\n (%d , %d, %d) , %d \n", width, height, color_type, bit_depth);
	//printf(" Size of %d \n", sizeof(png_byte));

	number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);


	/* read file */
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[read_png_file] Error during read_image");

	CShape s;
	s.width = width;
	s.height = height;
	
//	fprintf(stderr," Color Type = %d \n",(int)color_type);
	
	if (color_type == 0)
		s.nBands = 1;
	else if (color_type == 2)
		s.nBands = 3;  // setting for current RGB format png

	

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (y=0; y<height; y++)
		row_pointers[y] = (png_byte*) malloc(info_ptr->rowbytes);

	//printf(" Size of %d %d %d %d \n", sizeof(png_byte), sizeof(png_bytep)*height, sizeof(row_pointers), info_ptr->rowbytes);

	png_read_image(png_ptr, row_pointers);

	//printf("addr = %p %d", row_pointers, *row_pointers);

	img.ReAllocate(s, bit_depth/8, row_pointers);

  fclose(fp);

	

	/*
	if (info_ptr->color_type != PNG_COLOR_TYPE_RGB)
		abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (is %d)", info_ptr->color_type);

	for (y=0; y<height; y++) {
		png_byte* row = row_pointers[y];
		for (x=0; x<width; x++) {
			//png_byte* ptr = &(row[x*6]);
			unsigned short* ptr = (unsigned short*) &(row[x*6]);
			if (x==20 && y ==256)
				printf("Pixel at position [ %d - %d ] has the following RGB values: %d - %d - %d\n", x, y, endianb2l(ptr[0]), endianb2l(ptr[1]), endianb2l(ptr[2]));

		}
	}
	*/


}




void write_png_file(CImage2 img, char* file_name)
{

// color_type = 2 means 3 16 bit values for R, G and B respectively
// color_type = 0 means 1 16 bit value for grayscale

// for now i will set it to 0
	int color_type = 0, bit_depth = 16;

	/* create file */
	FILE *fp = fopen(file_name, "wb");
	if (!fp)
		abort_("[write_png_file] File %s could not be opened for writing", file_name);


	/* initialize stuff */
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	
	if (!png_ptr)
		abort_("[write_png_file] png_create_write_struct failed");

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
		abort_("[write_png_file] png_create_info_struct failed");

	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during init_io");

	png_init_io(png_ptr, fp);


	/* write header */
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during writing header");

	png_set_IHDR(png_ptr, info_ptr, img.Shape().width, img.Shape().height,
		     bit_depth, color_type, PNG_INTERLACE_NONE,
		     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);


	/* write bytes */
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during writing bytes");

	png_write_image(png_ptr, img.getmemStart());


	/* end write */
	if (setjmp(png_jmpbuf(png_ptr)))
		abort_("[write_png_file] Error during end of write");

	png_write_end(png_ptr, NULL);

        /* cleanup heap allocation */
	//for (y=0; y<height; y++)
	//	free(row_pointers[y]);
	//free(row_pointers);

        fclose(fp);
}





















void ReadFilePNG(CByteImage& img, const char* filename)
{
    // open the PNG input file
    FILE *stream = fopen(filename, "rb");
    if (stream == 0)
        throw CError("ReadFilePNG: could not open %s", filename);

    // first check the eight byte PNG signature
    png_byte pbSig[8];
    fread(pbSig, 1, 8, stream);
	if (!png_check_sig(pbSig, 8)) {
        fclose(stream);
        throw CError("ReadFilePNG: invalid PNG signature");
	}

    // create the two png(-info) structures
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
      (png_error_ptr)pngfile_error, (png_error_ptr)NULL);

	if (!png_ptr) {
        fclose(stream);
		throw CError("ReadFilePNG: error creating png structure");
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_read_struct(&png_ptr, NULL, NULL);
        fclose(stream);
		throw CError("ReadFilePNG: error creating png structure");
	}

	png_init_io(png_ptr, stream);
	png_set_sig_bytes(png_ptr, 8);

	// read all PNG info up to image data
	png_read_info(png_ptr, info_ptr);

	// get width, height, bit-depth and color-type
	int width, height, bits, colorType, nBands;

	png_get_IHDR(png_ptr, info_ptr, 
		(png_uint_32 *)&width, (png_uint_32 *)&height,
		&bits, &colorType, NULL, NULL, NULL);
	nBands = (int)png_get_channels(png_ptr, info_ptr);

	if (DEBUG_ImageIOpng)
	fprintf(stderr, " w=%d, h=%d, %2d bits, %s, nB=%d",
		width, height, bits,
		colorType == PNG_COLOR_TYPE_GRAY ? "gray" :
		colorType == PNG_COLOR_TYPE_PALETTE ? "plt " :
		colorType == PNG_COLOR_TYPE_RGB ? "rgb " :
		colorType == PNG_COLOR_TYPE_RGB_ALPHA ? "rgba" :
		colorType == PNG_COLOR_TYPE_GRAY_ALPHA ? "gr-a" : "??? ",
		nBands);


	// get rid of lower-order byte in 16-bit images
	// TODO: could allow this and read in IntImage in this case...
	if (bits == 16)
		png_set_strip_16(png_ptr);

	// change palette color into RGB
	if (colorType == PNG_COLOR_TYPE_PALETTE)
		png_set_expand(png_ptr);

	// want at least 8 bits
	if (bits < 8)
		png_set_expand(png_ptr);

	// if there is a transparent palette entry, create alpha channel
	if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
		png_set_expand(png_ptr);

	// make gray images with alpha channel into RGBA -- TODO: or just ignore alpha?
	if (colorType == PNG_COLOR_TYPE_GRAY_ALPHA)
		// colorType == PNG_COLOR_TYPE_GRAY       // but leave gray images alone
		png_set_gray_to_rgb(png_ptr);

	// set the background color to draw transparent and alpha images over.
	// only needed for gray images with alpha 
	if (colorType == PNG_COLOR_TYPE_GRAY_ALPHA ||
		colorType == PNG_COLOR_TYPE_GRAY) {
		png_color_16 *pBackground;
		if (png_get_bKGD(png_ptr, info_ptr, &pBackground))
			png_set_background(png_ptr, pBackground, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0);
		}

	// if required set gamma conversion
	// this seems to cause problems, so let's just leave gamma alone.
	//double gamma;
	//if (png_get_gAMA(png_ptr, info_ptr, &gamma)) {
	// //fprintf(stderr, "\n reading gamma %lf\n", gamma);
	//png_set_gamma(png_ptr, 1.0, gamma);
	//}

	// we need colors in BGR order, not RGB
	png_set_bgr(png_ptr);

	// always convert 3-band to 4-band images (add alpha):
	if (colorType == PNG_COLOR_TYPE_RGB)
		png_set_add_alpha(png_ptr, 255, PNG_FILLER_AFTER);

	// after the transformations have been registered update info_ptr data
	png_read_update_info(png_ptr, info_ptr);

	// get again width, height and the new bit-depth and color-type
	png_get_IHDR(png_ptr, info_ptr, 
		(png_uint_32 *)&width, (png_uint_32 *)&height,
		&bits, &colorType, NULL, NULL, NULL);
	nBands = (int)png_get_channels(png_ptr, info_ptr);

	if (DEBUG_ImageIOpng)
	fprintf(stderr, "  -> w=%d, h=%d, %2d bits, %s, nB=%d\n",
		width, height, bits,
		colorType == PNG_COLOR_TYPE_GRAY ? "gray" :
		colorType == PNG_COLOR_TYPE_PALETTE ? "plt " :
		colorType == PNG_COLOR_TYPE_RGB ? "rgb " :
		colorType == PNG_COLOR_TYPE_RGB_ALPHA ? "rgba" :
		colorType == PNG_COLOR_TYPE_GRAY_ALPHA ? "gr-a" : "??? ",
		nBands);
	

	if (! (nBands==1 || nBands==3 || nBands==4)) {
        fclose(stream);
		throw CError("ReadFilePNG: Can't handle nBands=%d", nBands);
	}

	// Set the image shape
	CShape sh(width, height, nBands);

	// Allocate the image if necessary
	img.ReAllocate(sh);

	//  allocate a vector of row pointers
	std::vector<uchar *> rowPtrs;
	rowPtrs.resize(height);
	for (int y = 0; y<height; y++)
		rowPtrs[y] = &img.Pixel(0, y, 0);

	// read the whole image
	png_read_image(png_ptr, &rowPtrs[0]);

 	// read the additional chunks in the PNG file (not really needed)
	png_read_end(png_ptr, NULL);

	png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

	fclose(stream);
}


void WriteFilePNG(CByteImage img, const char* filename)
{
	img = removeRedundantBands(img);

    CShape sh = img.Shape();
    int width = sh.width, height = sh.height, nBands = sh.nBands;

	// Make sure the image has the smallest number of bands before writing.
	// That is, if it's 4 bands with full alpha, reduce to 3 bands.  
	// If it's 3 bands with constant colors, make it 1-band.

    FILE *stream = fopen(filename, "wb");
    if (stream == 0)
        throw CError("WriteFilePNG: could not open %s", filename);

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
      (png_error_ptr)pngfile_error, (png_error_ptr)NULL);

	if (!png_ptr) {
        fclose(stream);
		throw CError("WriteFilePNG: error creating png structure");
	}

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        fclose(stream);
		throw CError("WriteFilePNG: error creating png structure");
    }

	png_init_io(png_ptr, stream);

	int bits = 8;
	int colortype =
		nBands == 1 ? PNG_COLOR_TYPE_GRAY :
		nBands == 3 ? PNG_COLOR_TYPE_RGB :
	                  PNG_COLOR_TYPE_RGB_ALPHA;
	png_set_IHDR(png_ptr, info_ptr, width, height, 
		bits, colortype,
		PNG_INTERLACE_NONE, 
		PNG_COMPRESSION_TYPE_DEFAULT, 
		PNG_FILTER_TYPE_DEFAULT);

	// write the file header information
	png_write_info(png_ptr, info_ptr);

	// swap the BGR pixels in the DiData structure to RGB
	png_set_bgr(png_ptr);

	//  allocate a vector of row pointers
	std::vector<uchar *> rowPtrs;
	rowPtrs.resize(height);
	for (int y = 0; y<height; y++)
		rowPtrs[y] = &img.Pixel(0, y, 0);

	// write the whole image
	png_write_image(png_ptr, &rowPtrs[0]);

	// write the additional chunks to the PNG file (not really needed)
	png_write_end(png_ptr, info_ptr);

	png_destroy_write_struct(&png_ptr, (png_infopp) NULL);

    fclose (stream);
}
