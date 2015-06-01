#include "sphere.h"

Sphere::Sphere(float* cent, float rad)
{
    this->center = new float[3];
    memcpy(center,cent,3*sizeof(float));
    this->radius = rad;
}
