/**
 * Pixel.h
 *
 * Pixel class.
 *
 * original by Wagner Correa at Princeton, 1999
 * turned to C++ by Robert Osada, 2000
 **/
#ifndef PIXEL_INCLUDED
#define PIXEL_INCLUDED


/**
 * Component, fundamental type
 **/
typedef unsigned char Component;

// Confines a component in the range [0..255]
inline Component ComponentClamp(int i)
{ return (i<0) ? 0 : (i>255) ? 255 : i; }

// Returns a random number in the range [0.255]
Component ComponentRandom(void);

// Scales the component by the given factor
Component ComponentScale(Component c, double f);

// Linear interpolation of the components
// Returns (1 - t) * c + t * d
Component ComponentLerp(Component c, Component d, double t);



/**
 * Pixel
 **/
struct Pixel
{
    // Data
    Component r, g, b, a;

    // Constructor
    Pixel (Component r_=0, Component g_=0, Component b_=0, Component a_=255) : r(r_), g(g_), b(b_), a(a_) {}

    // Set
    void Set (Component  r_, Component  g_, Component  b_, Component  a_) { r=r_; g=g_; b=b_; a=a_; }
    void Set (Component  r_, Component  g_, Component  b_)                { r=r_; g=g_; b=b_; }

    void SetClamp (double r_, double g_, double b_);
    void SetClamp (double r_, double g_, double b_, double a_);

    // Returns the luminance of the pixel.
    Component Luminance ();
};

// Returns a random pixel.
Pixel PixelRandom (void);

// Component-wise addition of pixels.
Pixel operator+ (const Pixel& p, const Pixel& q);

// Component-wise multiplication of pixels.
Pixel operator* (const Pixel& p, const Pixel& q);

// Component-wise multiplication of pixel by scalar.
Pixel operator* (const Pixel& p, double f);

// Linear interpolation of pixels.  Returns (1 - t) * p + t * q.
Pixel PixelLerp (const Pixel& p, const Pixel& q, double t);

/**
 * Composites two pixels.  Returns beta * top + (1 - beta) * alpha *
 * bottom, where alpha and beta are the alpha components of bottom and
 * top scaled to the [0..1] range.
 **/
Pixel PixelOver (const Pixel& bottom, const Pixel& top);


#endif

