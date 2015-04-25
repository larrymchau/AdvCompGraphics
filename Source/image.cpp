#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <time.h>


#define M_PI 3.14159265358979323846
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

void Image::Blur(int n)
{
  /* Your Work Here (Section 3.4.1) */
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			float r = 0; float g = 0; float b = 0;
			float weightSum = 0;

			for (int tempX = i - (n/2); tempX <= i + (n/2); tempX++){
				for (int tempY = j - (n/2); tempY <= j + (n/2); tempY++){
					if (tempX < 0 || tempX > width - 1){
						continue;
					}
					if (tempY < 0 || tempY > height - 1){
						continue;
					}
					float sigma = floorf(n / 2) / 2;
					float weight = exp(-((tempX - i)*(tempX - i) + (tempY - j)*(tempY - j))/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
					//fprintf(stderr,"%d %d\n", tempX , tempY);
					r += pixels[tempX*height + tempY].r * weight;
					g += pixels[tempX*height + tempY].g * weight;
					b += pixels[tempX*height + tempY].b * weight;
					weightSum += weight;
				}
			}
			r = (int)((r / weightSum) + 0.5);
			g = (int)((g / weightSum) + 0.5);
			b = (int)((b / weightSum) + 0.5);
			pixels[i*height + j].SetClamp(r, g, b);
		}
	}

}

void Image::Sharpen() 
{
	//int normalizeFactor = 7;
	//std::cout << "into sharpen" << std::endl;
	//printf("%s\n", "into sharpen");
	float* sharpenWeight = new float[9];
	sharpenWeight[0] = -1/7.0f; sharpenWeight[1] = -2/7.0f; sharpenWeight[2] = -1/7.0f;
	sharpenWeight[3] = -2/7.0f; sharpenWeight[4] = 19/7.0f; sharpenWeight[5] = -2/7.0f;
	sharpenWeight[6] = -1/7.0f; sharpenWeight[7] = -2/7.0f; sharpenWeight[8] = -1/7.0f;
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			float r = 0; float g = 0; float b = 0;
			float weightSum = 0;
			int count = 0;
			for (int tempX = i - 1; tempX <= i + 1; tempX++){
				for (int tempY = j - 1; tempY <= j + 1; tempY++,count++){
					if (tempX < 0 || tempX > width - 1){
						continue;
					}
					if (tempY < 0 || tempY > height - 1){
						continue;
					}
					float weight = sharpenWeight[count];
					r += pixels[tempX*height + tempY].r * weight;
					g += pixels[tempX*height + tempY].g * weight;
					b += pixels[tempX*height + tempY].b * weight;
					weightSum += weight;
				}
			}
			r = (int)((r / weightSum) + 0.5);
			g = (int)((g / weightSum) + 0.5);
			b = (int)((b / weightSum) + 0.5);
			pixels[i*height + j].SetClamp(r, g, b);
		}
	}

  /* Your Work Here (Section 3.4.2) */
}

void Image::EdgeDetect(int threshold)
{
  /* Your Work Here (Section 3.4.3) */
	//std::cout << "into edgedetect" << std::endl;
	float* horizontalWeight = new float[9];
	horizontalWeight[0] = -1; horizontalWeight[1] = 0; horizontalWeight[2] = 1;
	horizontalWeight[3] = -2; horizontalWeight[4] = 0; horizontalWeight[5] = 2;
	horizontalWeight[6] = -1; horizontalWeight[7] = 0; horizontalWeight[8] = 1;

	float* verticalWeight = new float[9];
	verticalWeight[0] = 1; verticalWeight[1] = 2; verticalWeight[2] = 1;
	verticalWeight[3] = 0; verticalWeight[4] = 0; verticalWeight[5] = 0;
	verticalWeight[6] = -1; verticalWeight[7] = -2; verticalWeight[8] = -1;

	//Pixel *hPixels = new Pixel[num_pixels];
	//Pixel *vPixels = new Pixel[num_pixels];
	float *hPixelsr = new float[num_pixels];
	float *hPixelsg = new float[num_pixels];
	float *hPixelsb = new float[num_pixels];

	float *vPixelsr = new float[num_pixels];
	float *vPixelsg = new float[num_pixels];
	float *vPixelsb = new float[num_pixels];

	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			float hr = 0; float hg = 0; float hb = 0;
			//float hweightSum = 0;
			float vr = 0; float vg = 0; float vb = 0;
			//float vweightSum = 0;
			int count = 0;
			//std::cout << "start loop for 1 pixel" << std::endl;
			for (int tempY = j - 1; tempY <= j + 1; tempY++){
				for (int tempX = i - 1; tempX <= i + 1; tempX++, count++){
					if (tempX < 0 || tempX > width - 1){
						continue;
					}
					if (tempY < 0 || tempY > height - 1){
						continue;
					}
					float weight = horizontalWeight[count];
					//fprintf(stderr, "horizontal count: %d : %.2f // tempX: %d , tempY: %d, x: %d, y: %d\n", count, weight,tempX,tempY,i,j);
					hr += pixels[tempX*height + tempY].r * weight;
					hg += pixels[tempX*height + tempY].g * weight;
					hb += pixels[tempX*height + tempY].b * weight;
					//hweightSum += weight;
					
					weight = verticalWeight[count];
					//fprintf(stderr, "vertical count: %d : %.2f // tempX: %d , tempY: %d, x: %d, y: %d\n", count, weight, tempX, tempY, i, j);
					vr += pixels[tempX*height + tempY].r * weight;
					vg += pixels[tempX*height + tempY].g * weight;
					vb += pixels[tempX*height + tempY].b * weight;
					//vweightSum += weight;
				}
			}
			
			//std::cout << "end loop for 1 pixel" << std::endl;
			hPixelsr[i*height + j] = (hr);// / hweightSum);
			hPixelsg[i*height + j] = (hg);// / hweightSum);
			hPixelsb[i*height + j] = (hb);// / hweightSum);
			
			vPixelsr[i*height + j] = (vr);// / vweightSum);
			vPixelsg[i*height + j] = (vg);// / vweightSum);
			vPixelsb[i*height + j] = (vb);// / vweightSum);
			//r * 76 + g * 150 + b * 29
			//Pixel h = Pixel(hr, hg, hb);
			//Pixel v = Pixel(vr, vg, vb);
			float hPixelLuminance = ((76.0f / 255.0)*hPixelsr[i*height + j] + (150.0f / 255.0f)*hPixelsg[i*height + j] + (29.0f / 255.0f)*hPixelsb[i*height + j]);
			float vPixelLuminance = ((76.0f / 255.0)*vPixelsr[i*height + j] + (150.0f / 255.0f)*vPixelsg[i*height + j] + (29.0f / 255.0f)*vPixelsb[i*height + j]);
			float gradient = sqrtf((pow(abs(hPixelLuminance), 2)) + (pow(abs(vPixelLuminance), 2)));
			//float gradient = sqrt((h.Luminance() * h.Luminance()) + (v.Luminance() + v.Luminance()));
			//fprintf(stderr, "%.6f %.6f %.6f\n", hPixelsr[i*height + j], hPixelsg[i*height + j], hPixelsb[i*height + j]);
			//fprintf(stderr, "%.6f\n", gradient);
			//	fprintf(stderr, "%.6f %.6f %.6f %.6f\n", hr, hg, hb, hweightSum);
			//float gradient = sqrtf((hPixels[i*width + j].Luminance() * hPixels[i*width + j].Luminance()) + (vPixels[i*width + j].Luminance() * vPixels[i*width + j].Luminance()));
			if (gradient >= threshold)
				pixels[i*height + j].SetClamp(255,255,255);
			else{
				pixels[i*height + j].SetClamp(0, 0, 0);
			}
		}
	}
}


Image* Image::Scale(int sizex, int sizey)
{
  /* Your Work Here (Section 3.5.1) */
  double sx = (double(sizex)) / (double(width));
  double sy = (double(sizey)) / (double(height));

  Image* img = new Image(sizex, sizey);
  Pixel* pix = img->pixels;

  Pixel* newPixels = new Pixel[sizex * sizey];
	
  for (int i = 0; i < sizex; i++) {
	  for (int j = 0; j < sizey; j++){
		  newPixels[(j*sizex) + i] = Sample(i, j, sx, sy, 0);
	  }
  }

  //pix[(j*sizex) + i] = Sample(i, j, sx, sy, 0, 1);
  width = sizex;
  height = sizey;
  num_pixels = sizex * sizey;
  pixels = newPixels;

  //Pixel* newPixels2 = new Pixel[sizex * sizey];
  

  for (int i = 0; i < sizex; i++) {
	  for (int j = 0; j < sizey; j++){
		  pix[(j*sizex) + i] = Sample(i, j, sx, sy, 0, 1);
	  }
  }

  return  img;
}

void Image::Shift(double sx, double sy)
{
	Pixel* newPixels = new Pixel[num_pixels];
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++){
			newPixels[(j*width) + i] = Sample(i, j, sx, sy, 1);
		}
	}
	delete pixels;
	pixels = newPixels;

	return;
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

Pixel Image::Sample (double u, double v, double sx, double sy, bool shift, int pass)
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
	  if (shift){
		  int px = ceil(u - sx);
		  int py = ceil(v - sy);
		  if (px < 0 || px > width) {
			  return Pixel();
		  }
		  if (py < 0 || py > height){
			  return Pixel();
		  }
		  return pixels[py*width + px];
	  }
	  else{
		  int px = u / sx;
		  int py = v / sy;
		  return pixels[py*width + px];
	  }
  }

  else if (sampling_method == IMAGE_SAMPLING_HAT) {
	  if (shift){
		  return ConvolveHat(1, 1, u, v, sx, sy, 1);
	  }
	  else {
		  if (sx <= 1 && sy <= 1){
			  return ConvolveHat(sx, sy, u, v, sx, sy, 0);
		  }
		  if (sx <= 1 && sy >= 1){
			  return ConvolveHat(sx, 1, u, v, sx, sy, 0);
		  }
		  if (sx >= 1 && sy <= 1){
			  return ConvolveHat(1, sy, u, v, sx, sy, 0);
		  }
		  if (sx >= 1 && sy >= 1) {
			  return ConvolveHat(1, 1, u, v, sx, sy, 0);
		  }
		  return pixels[(int)(v*width + u)];
	  }
  }

  else if (sampling_method == IMAGE_SAMPLING_MITCHELL) {
	  if (shift) {
		  return ConvolveMitchell(1, 1, u, v, sx, sy, 1,pass);
	  }
	  else {
		  if (sx <= 1 && sy <= 1){
			  return ConvolveMitchell(sx, sy, u, v, sx, sy, 0,pass);
		  }
		  if (sx <= 1 && sy >= 1){
			  return ConvolveMitchell(sx, 1, u, v, sx, sy, 0,pass);
		  }
		  if (sx >= 1 && sy <= 1){
			  return ConvolveMitchell(1, sy, u, v, sx, sy, 0,pass);
		  }
		  if (sx >= 1 && sy >= 1) {
			  return ConvolveMitchell(1, 1, u, v, sx, sy, 0,pass);
		  }
		  return pixels[(int)(v*width + u)];
	  }
  }

  else {
    fprintf(stderr,"I don't understand what sampling method is used\n") ;
    exit(1) ;
  }

  return Pixel();
}

Pixel Image::ConvolveHat(double paramX, double paramY, double u, double v, double sx, double sy, bool shift) {
	double px, py;
	double red = 0; double green = 0; double blue = 0;
	double normalize = 0;
	bool Xdir = true; bool Ydir = true;
	if (shift){
		
		if (sx >= 0) //shifting right
			Xdir = true;
		else  //shifting left
			Xdir = false;
		if (sy >= 0)
			Ydir = true;
		else
			Ydir = false;
		for (int x = ceil(u - sx - 1); (x <= u - sx + 1);  x++ ){
			for (int y = ceil(v - sy - 1); (y <= v - sy + 1) ; y++ ){
				
				if (Xdir){
					if (x < 0) {
						//return Pixel();//continue;
						continue;
					}
					else if(x >= width){
						continue;
					}
				}
				else{
					if (x >= width){
						//return Pixel();
						continue;
					}
					else if (x < 0){
						continue;
					}
				}
				
				if (Ydir){
					if (y < 0) {
						//return Pixel();//continue;
						continue;
					}
					else if (y >= height){
						continue;
					}
				}
				else{
					if (y >= height){
						//return Pixel();
						continue;
					}
					else if (y < 0){
						continue;
					}
				}
				//std::cerr << "x: " << x << " y: " << y << std::endl;
				px = Barlett(((double)x) - u +sx);
				py = Barlett(((double)y) - v +sy);
				red += (double)pixels[(int)(y*width + x)].r * px * py;
				green += (double)pixels[(int)(y*width + x)].g * px * py;
				blue += (double)pixels[(int)(y*width + x)].b * px * py;
				normalize += px * py;
			}
		}
		return Pixel(red / normalize, green / normalize, blue / normalize);
	}
	else{
		for (int x = (u - 1) / sx; x <= (u + 1) / sx; x++){
			for (int y = (v - 1) / sy; y <= (v + 1) / sy; y++) {
				if (x < 0 || x >= width || y < 0 || y >= height) {
					continue;
				}
				px = Barlett(paramX*(((double)x) - (u / sx)));
				py = Barlett(paramY*(((double)y) - (v / sy)));
				red += (double)pixels[(int)(y*width + x)].r * px * py;
				green += (double)pixels[(int)(y*width + x)].g * px * py;
				blue += (double)pixels[(int)(y*width + x)].b * px * py;
				normalize += px * py;
			}
		}

		return Pixel(red / normalize, green / normalize, blue / normalize);
	}
}

Pixel Image::ConvolveMitchell(double paramX, double paramY, double u, double v, double sx, double sy, bool shift, int pass) {
	double px, py;
	double red = 0; double green = 0; double blue = 0;
	double normalize = 0;
	/*
		
	*/
	bool Xdir = true; bool Ydir = true;
	if (shift){
		if (sx >= 0) //shifting right
			Xdir = true;
		else  //shifting left
			Xdir = false;
		if (sy >= 0)
			Ydir = true;
		else
			Ydir = false;
		
		for (int x = (u - sx - 2); (x < u - sx + 1); x++){
			for (int y = (v - sy - 2); (y < v - sy + 1); y++){

				
				if (Xdir){
					if (x < 0) {
						//return Pixel();//continue;
						continue;
					}
					else if (x >= width){
						continue;
					}
				}
				else{
					if (x >= width){
						return Pixel();
						continue;
					}
					else if (x < 0){
						continue;
					}
				}

				if (Ydir){
					if (y < 0) {
						//return Pixel();//continue;
						continue;
					}
					else if (y >= height){
						continue;
					}
				}
				else{
					if (y >= height){
						return Pixel();
						continue;
					}
					else if (y < 0){
						continue;
					}
				}
				//std::cerr << "x: " << x << " y: " << y << std::endl;
				px = Mitchell((((double)x) - u +sx));
				py = Mitchell((((double)y) - v +sy));
				red += (double)pixels[(int)(y*width + x)].r * px * py;
				green += (double)pixels[(int)(y*width + x)].g * px * py;
				blue += (double)pixels[(int)(y*width + x)].b * px * py;
				normalize += px * py;
			}
		}

		return Pixel(red / normalize, green / normalize, blue / normalize);

	}
	
	else{
		int count = 0;
		for (int x = ceil(((double)u / sx) - 2); x < (((double)u / sx) + 2); x++,count=0){
			for (int y = ceil(((double)v / sy) - 2); y < (((double)v / sy) + 2); y++,count++) {
				//if (count > 3){
					//std::cerr << count << std::endl;
				//}

				//std::cerr << "x: " << x << " y: " << y <<" "; 
				if (x < 0 || x >= width || y < 0 || y >= height) {
					continue;
				}
				px = Mitchell(paramX*(((double)x) - ((double)u / sx)));
				py = Mitchell(paramY*(((double)y) - ((double)v / sy)));
				if (pass == 0){
					red += (double)pixels[(int)(y*width + x)].r * px;
					green += (double)pixels[(int)(y*width + x)].g * px;
					blue += (double)pixels[(int)(y*width + x)].b * px;
					normalize += px;
				}
				else{
					red += (double)pixels[(int)(y*width + x)].r * py;
					green += (double)pixels[(int)(y*width + x)].g * py;
					blue += (double)pixels[(int)(y*width + x)].b * py;
					normalize += py;
				}
				//if (x == 1 && y == 6)
					//std::cerr << "px: " << px << " py: " << py << " normalize: " << normalize << std::endl;
			}
		}

		//std::cerr << std::endl;


		return Pixel(red / normalize, green / normalize, blue / normalize);
	}
}

double Image::Barlett(double x){
	double temp = abs(x);
	if (temp > 1){
		return 0;
	}
	return 1 - temp;
}

double Image::Mitchell(double x){
	double temp = abs(x);
	if (0 <= temp && temp < 1)
	{
		return (1.0 / 6.0)*(7.0*pow(temp, 3) - 12.0*pow(temp, 2) + 16.0 / 3.0);
	}
	else if (1 <= temp && temp < 2)
	{
		return (1.0 / 6.0)*((-7.0 / 3.0)*pow(temp, 3) + 12 * pow(temp, 2) - 20 * temp + 32.0 / 3.0);
	}
	return 0;
}
