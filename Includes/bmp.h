/**
 * bmp.h
 *
 * reads and writes 24-bit uncompressed bitmaps
 **/
#ifndef BMP_INCLUDED
#define BMP_INCLUDED

#include "image.h"


// Reads an image from a BMP file.
Image *BMPReadImage(FILE *fp);

// Writes an image to a BMP file.
void BMPWriteImage(Image *img, FILE *fp);


#endif
