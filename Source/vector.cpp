#include "vector.h"

Vector Vector::Rotate (double theta) const
{
    double c = cos(theta);
    double s = sin(theta);

    return Vector(x*c - y*s, x*s + y*c);
}

