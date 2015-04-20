#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <time.h>

/**
 * Image
 **/
Image::Image (int width_, int height_)
{
    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
}


Image::Image (const Image& src)
{
    width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
    memcpy(pixels, src.pixels, src.width * src.height * sizeof(Pixel));
}


Image::~Image ()
{
    delete [] pixels;
    pixels = NULL;
}

/*
void Image::AddNoise (double factor)
{

}
*/

void Image::Brighten(double factor)
{
	/* Your Work Here  (section 3.2.1 of assignment)*/
	if (factor < 0){
		factor = 0.0;
	}
	for (int i = 0; i < num_pixels; i++){
		pixels[i].SetClamp(pixels[i].r * factor, pixels[i].g * factor, pixels[i].b * factor);
	}
}


void Image::ChangeContrast(double factor)
{
	/* Your Work Here (section 3.2.2) */
	double avg = 0.0;
	for (int i = 0; i < num_pixels; i++){
		avg += pixels[i].Luminance();
	}
	avg /= num_pixels;
	Pixel grey;
	grey.SetClamp(avg, avg, avg);
	for (int i = 0; i < num_pixels; i++){
		pixels[i] = PixelLerp(grey, pixels[i], factor);
	}
}


void Image::ChangeSaturation(double factor)
{
	/* Your Work Here (section 3.2.3) */
	for (int i = 0; i < num_pixels; i++){
		double greyVal = pixels[i].Luminance();
		Pixel grey; grey.SetClamp(greyVal, greyVal, greyVal);
		pixels[i] = PixelLerp(grey, pixels[i], factor);
	}
}

void Image::ChangeGamma(double factor)
{
	/* Your Work Here (section 3.2.4) */
	double red, green, blue;
	for (int i = 0; i < num_pixels; i++){
		//convert pixel [0-255] to [0-1] -> Cold
		//then do Cnew = Cold^(1/y) and then multiply by 255 again to get
		//[0-255] domain
		red = pow((double)(pixels[i].r/255.0), ((double)1.0) / factor)*255.0;
		green = pow((double)(pixels[i].g/255.0), ((double)1.0) / factor)*255.0;
		blue = pow((double)(pixels[i].b/255.0), ((double)1.0) / factor)*255.0;
		this->pixels[i].SetClamp(red,green,blue);
    }
}

Image* Image::Crop(int x, int y, int w, int h)
{
  /* Your Work Here (section 3.2.5) */
	if (x + w > width || y + w > height){
		std::cerr << "Out of Bounds" << std::endl;
		return NULL;
	}
	if (x < 0 || y < 0 || w < 0 || h < 0) {
		std::cerr << "Invalid Parameters" << std::endl;
		return NULL;
	}
	Image* img = new Image(w, h);
	Pixel* croppedImage = img->pixels;
	int cIdx = 0;
	//start position
	int Idx = y*(width)+x;
	/* row major */
	for (int i = 0; i < h; i = i ++){
		for (int j = 0; j < w; j++){
			croppedImage[cIdx] = pixels[Idx + j];
			cIdx++;
		}
		Idx += width;
	}

	return img;
}

/*
void Image::ExtractChannel(int channel)
{
  // For extracting a channel (R,G,B) of image.  
  // Not required for the assignment
}
*/

void Image::Quantize (int nbits)
{
	double b = pow(2, nbits);
  /* Your Work Here (Section 3.3.1) */
	double red, green, blue;
	for (int i = 0; i < num_pixels; i++) {
		double factor = pow(2, nbits);
		//cf = floor(255*floor(p*b)/(b-1))
		red = floor(255.0*floor(pixels[i].r*factor / 256.0) / (factor - 1));
		green = floor(255.0*floor(pixels[i].g*factor / 256.0) / (factor - 1));
		blue = floor(255.0*floor(pixels[i].b*factor / 256.0) / (factor - 1));
		pixels[i].SetClamp(red, green, blue);
	}
}


void Image::RandomDither (int nbits)
{
  /* Your Work Here (Section 3.3.2) */

	srand(time(NULL));
	double red, green, blue, random;
	for (int i = 0; i < num_pixels; i++) {
		double factor = pow(2, nbits);
		//cf = floor(255*floor(p*b)/(b-1))
		random = ((double)rand() / (double)RAND_MAX) - 0.5; //rand number -0.5-0.5
		red = floor(255.0*floor((pixels[i].r*factor / 256.0) + random) / (factor - 1));
		green = floor(255.0*floor((pixels[i].g*factor / 256.0) + random) / (factor - 1));
		blue = floor(255.0*floor((pixels[i].b*factor / 256.0) + random) / (factor - 1));
		pixels[i].SetClamp(red, green, blue);
	}
}


/* Matrix for Bayer's 4x4 pattern dither. */
/* uncomment its definition if you need it */

/*
static int Bayer4[4][4] =
{
    {15, 7, 13, 5},
    {3, 11, 1, 9},
    {12, 4, 14, 6},
    {0, 8, 2, 10}
};


void Image::OrderedDither(int nbits)
{
  // For ordered dithering
  // Not required for the assignment
}

*/

/* Error-diffusion parameters for Floyd-Steinberg*/
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	/* Your Work Here (Section 3.3.3) */

	double red, green, blue;
	for (int i = 0; i < num_pixels; i++) {
		double factor = pow(2, nbits);
		//cf = floor(255*floor(p*b)/(b-1))
		red = floor(255.0*floor(pixels[i].r*factor / 256.0) / (factor - 1));
		green = floor(255.0*floor(pixels[i].g*factor / 256.0) / (factor - 1));
		blue = floor(255.0*floor(pixels[i].b*factor / 256.0) / (factor - 1));
		if (i < num_pixels - 1){
			pixels[i+1].r += (7.0 / 16.0) * (pixels[i].r - red);
			pixels[i+1].g += (7.0 / 16.0) * (pixels[i].g - green);
			pixels[i+1].b += (7.0 / 16.0) * (pixels[i].b - blue);
		}
		if (i < num_pixels - width){
			pixels[i+width-1].r += (3.0 / 16.0) * (pixels[i].r - red);
			pixels[i+width-1].g += (3.0 / 16.0) * (pixels[i].g - green);
			pixels[i+width-1].b += (3.0 / 16.0) * (pixels[i].b - blue);

			pixels[i+width].r += (5.0 / 16.0) * (pixels[i].r - red);
			pixels[i+width].g += (5.0 / 16.0) * (pixels[i].g - green);
			pixels[i+width].b += (5.0 / 16.0) * (pixels[i].b - blue);
		}
		if (i < num_pixels - width - 1){
			pixels[i + width - 1].r += (1.0 / 16.0) * (pixels[i].r - red);
			pixels[i + width - 1].g += (1.0 / 16.0) * (pixels[i].g - green);
			pixels[i + width - 1].b += (1.0 / 16.0) * (pixels[i].b - blue);
		}
		pixels[i].SetClamp(red, green, blue);
	}

}

void ImageComposite(Image *bottom, Image *top, Image *result)
{
  // Extra Credit (Section 3.7).
  // This hook just takes the top image and bottom image, producing a result
  // You might want to define a series of compositing modes as OpenGL does
  // You will have to use the alpha channel here to create Mattes
  // One idea is to composite your face into a famous picture
}

void Image::Convolve(int *filter, int n, int normalization, int absval) {
  // This is my definition of an auxiliary function for image convolution 
  // with an integer filter of width n and certain normalization.
  // The absval param is to consider absolute values for edge detection.
  
  // It is helpful if you write an auxiliary convolve function.
  // But this form is just for guidance and is completely optional.
  // Your solution NEED NOT fill in this function at all
  // Or it can use an alternate form or definition
}

void Image::Blur(int n)
{
  /* Your Work Here (Section 3.4.1) */
}

void Image::Sharpen() 
{
  /* Your Work Here (Section 3.4.2) */
}

void Image::EdgeDetect(int threshold)
{
  /* Your Work Here (Section 3.4.3) */
}


Image* Image::Scale(int sizex, int sizey)
{
  /* Your Work Here (Section 3.5.1) */
  return NULL ;
}

void Image::Shift(double sx, double sy)
{
  /* Your Work Here (Section 3.5.2) */
}


/*
Image* Image::Rotate(double angle)
{
  // For rotation of the image
  // Not required in the assignment
  // But you can earn limited extra credit if you fill it in
  // (It isn't really that hard) 

    return NULL;
}
*/


void Image::Fun()
{
    /* Your Work Here (Section 3.6) */
}


Image* ImageMorph (Image* I0, Image* I1, int numLines, Line* L0, Line* L1, double t)
{
  /* Your Work Here (Section 3.7) */
  // This is extra credit.
  // You can modify the function definition. 
  // This definition takes two images I0 and I1, the number of lines for 
  // morphing, and a definition of corresponding line segments L0 and L1
  // t is a parameter ranging from 0 to 1.
  // For full credit, you must write a user interface to join corresponding 
  // lines.
  // As well as prepare movies 
  // An interactive slider to look at various morph positions would be good.
  // From Beier-Neely's SIGGRAPH 92 paper

    return NULL;
}


/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
  // Sets the filter to use for Scale and Shift
  // You need to implement point sampling, hat filter and mitchell

    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}

Pixel Image::Sample (double u, double v, double sx, double sy)
{
  // To sample the image in scale and shift
  // This is an auxiliary function that it is not essential you fill in or 
  // you may define it differently.
  // u and v are the floating point coords of the points to be sampled.
  // sx and sy correspond to the scale values. 
  // In the assignment, it says implement MinifyX MinifyY MagnifyX MagnifyY
  // separately.  That may be a better way to do it.
  // This hook is primarily to get you thinking about that you have to have 
  // some equivalent of this function.

  if (sampling_method == IMAGE_SAMPLING_POINT) {
    // Your work here
  }

  else if (sampling_method == IMAGE_SAMPLING_HAT) {
    // Your work here
  }

  else if (sampling_method == IMAGE_SAMPLING_MITCHELL) {
    // Your work here
  }

  else {
    fprintf(stderr,"I don't understand what sampling method is used\n") ;
    exit(1) ;
  }

  return Pixel() ;
}

