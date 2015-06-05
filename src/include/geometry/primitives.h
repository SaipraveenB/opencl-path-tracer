#include <stdgl.h>

#ifndef PRIM_H
#define PRIM_H

class Sphere
{
public:
  glm::vec4 center;
  float radius;

  Sphere(glm::vec4 cent, float rad);
  Sphere();
};

class Triangle
{
public:
  glm::vec4 p0;
  glm::vec4 p1;
  glm::vec4 p2;

  Triangle(glm::vec4 p0, glm::vec4 p1, glm::vec4 p2);
  Triangle();
};

class Plane
{
public:
  glm::vec4 normal;
  glm::vec4 point;

  Plane(glm::vec4 normal, glm::vec4 point);
  Plane();
};

class GeometryDescriptor{
public:
  int numSpheres;
  int numTriangles;
  int numPlanes;

  GeometryDescriptor( int a, int b, int c );
};


#endif
