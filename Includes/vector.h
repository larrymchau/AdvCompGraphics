#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

#include <math.h>

struct Vector
{
    double x, y;

    // constructor
    Vector (double x_=0, double y_=0) { x=x_; y=y_; }

    // set
    void Set (double x_, double y_) { x=x_; y=y_; }

    // length to origin
    double Len () const { return sqrt(x*x+y*y); }
    double LenSquare () const { return x*x+y*y; }

    // rotate vector around origin
    Vector Rotate (double theta) const;

    // return vector perpendicular to this one
    Vector Perp () const { return Vector(-y, x); }
};


// basic vector operations
inline Vector operator+ (const Vector& lhs, const Vector& rhs)
{
    return Vector(lhs.x+rhs.x, lhs.y+rhs.y);
}

inline Vector operator- (const Vector& lhs, const Vector& rhs)
{
    return Vector(lhs.x-rhs.x, lhs.y-rhs.y);
}

inline Vector operator* (const Vector& lhs, double rhs)
{
    return Vector(lhs.x*rhs, lhs.y*rhs);
}


// dot product
inline double Dot (const Vector& lhs, const Vector& rhs)
{
    return lhs.x*rhs.x + lhs.y*rhs.y;
}



#endif
