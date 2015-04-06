/**
 * Image.h
 *
 * Image class.
 *
 * original by Wagner Correa at Princeton, 1999
 * turned to C++ by Robert Osada, 2000
 * Slight modifications for 4162 assignment in Spring 2005.
 **/
#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED

#include <assert.h>
#include <stdio.h>
#include "pixel.h"
#include "vector.h"


/**
 * constants
 **/
enum {
    IMAGE_SAMPLING_POINT,
    IMAGE_SAMPLING_HAT,
    IMAGE_SAMPLING_MITCHELL,
    IMAGE_N_SAMPLING_METHODS
};

enum {
    IMAGE_CHANNEL_RED,
    IMAGE_CHANNEL_GREEN,
    IMAGE_CHANNEL_BLUE,
    IMAGE_CHANNEL_ALPHA,
    IMAGE_N_CHANNELS
};


/**
 * Line defined as two points (vectors).
 **/
struct Line
{
    Vector p, q;
};



/**
 * Image
 **/
class Image
{
public:
    Pixel *pixels;
    int width, height, num_pixels;
    int sampling_method;

public:
    // Creates a blank image with the given dimensions
    Image (int width, int height);

    // Copy iamage
    Image (const Image& src);

    // Destructor
    ~Image ();

    // Pixel access
    int ValidCoord (int x, int y)  { return x>=0 && x<width && y>=0 && y<height; }
    Pixel& GetPixel (int x, int y) { assert(ValidCoord(x,y)); return pixels[y*width + x]; }

    // Dimension access
    int Width     () { return width; }
    int Height    () { return height; }
    int NumPixels () { return num_pixels; }

    /*
    // Adds noise to an image.  The amount of noise is given by the factor
    // in the range [0.0..1.0].  0.0 adds no noise.  1.0 adds a lot of noise.
    void AddNoise (double factor);
    */

    // Brightens the image by multiplying each pixel component by the factor.
    void Brighten (double factor);

    /**
     * Changes the contrast of an image by interpolating between the image
     * and a constant gray image with the average luminance.
     * Interpolation reduces constrast, extrapolation boosts constrast,
     * and negative factors generate inverted images.
     **/
    void ChangeContrast (double factor);

    /**
     * Changes the saturation of an image by interpolating between the
     * image and a gray level version of the image.  Interpolation
     * decreases saturation, extrapolation increases it, negative factors
     * presrve luminance but invert the hue of the input image.
     **/
    void ChangeSaturation (double factor);

    /**
     * Gamma corrects the image, as defined for standard gamma correction
     * That is, p = p^(1/gamma) where p is the floating point color val
     **/
    void ChangeGamma (double factor) ;

    /**
     * Extracts a sub image from the image, at position (x, y), width w,
     * and height h.
     **/
    Image* Crop(int x, int y, int w, int h);

    /**
     * Extracts a channel of an image.  Leaves the specified channel
     * intact.  Sets all other ones to zero.
     **/
    //void ExtractChannel(int channel);

    /**
     * Quantizes an image with "nbits" bits per channel.
     **/
    void Quantize (int nbits);

    // Converts and image to nbits per channel using random dither.
    void RandomDither(int nbits);

    // Helper function to implement convolution 
    void Convolve(int *filter, int n, int normalization, int absval = 0) ;

    // Blurs an image with an n x n Gaussian filter.
    void Blur(int n);

    // Sharpens an image.
    void Sharpen(void);

    // Detects edges in an image.
    void EdgeDetect(int threshold);

    /**
     * Converts an image to nbits per channel using ordered dither, with a
     * 4x4 Bayer's pattern matrix.
     **/
    void OrderedDither(int nbits);

    /**
     * Converts an image to nbits per channel using Floyd-Steinberg dither
     * with error diffusion.
     **/
    void FloydSteinbergDither(int nbits);

    // Scales an image in x by sx, and y by sy.
    Image* Scale(int sizex, int sizey);

    // Translates an image in x by sx, and y by sy.
    void Shift(double sx, double sy);

    // Rotates an image by the given angle.
    //Image* Rotate(double angle);

    // Warps an image using a creative filter of your choice.
    void Fun();

    // Sets the sampling method.
    void SetSamplingMethod(int method);

    // Sample image using current sampling method.
    Pixel Sample(double u, double v, double sx, double sy);
};

// Composites the bottom and top images into the result image.
void ImageComposite(Image *bottom, Image *top, Image *result);

/**
 * Morph two images using [Beier92].
 * I0 and I1 are before and after images.
 * L0 and L1 are corresponding line segments.
 * t is the morth time, between 0 and 1.
 **/
Image* ImageMorph (Image* I0, Image* I1, int numLines, Line* L0, Line* L1, double t);


#endif
