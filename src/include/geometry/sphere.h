#include <stdgl.h>

#ifndef SPHERE_H
#define SPHERE_H

class Sphere
{
public:
    glm::vec4 center;
    float radius;

    Sphere(glm::vec4 cent, float rad);
    Sphere();
};


#endif
