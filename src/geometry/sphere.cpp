#include <geometry/sphere.h>


Sphere::Sphere(glm::vec4 center, float rad)
{
    this->center = center;
    this->radius = rad;
}

Sphere::Sphere(){
    this->center = glm::vec4(0.0f,0.0f,0.0f,0.0f);
    this->radius = 1.0f;
}
