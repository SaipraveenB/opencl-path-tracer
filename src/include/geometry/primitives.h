#include <stdgl.h>

#ifndef PRIM_H
#define PRIM_H

class Sphere
{
public:
	glm::vec4 center;
	float radius;
	unsigned int uSurf;
	unsigned int reserved0;
	unsigned int reserved1;
                    
	Sphere(glm::vec4 cent, float rad);
	Sphere();
};

class Triangle
{
public:
	glm::vec4 p0;
	glm::vec4 p1;
	glm::vec4 p2;
	unsigned int uSurf;
	unsigned int uReserved0;
	unsigned int uReserved1;
	unsigned int uReserved2;

	Triangle(glm::vec4 p0, glm::vec4 p1, glm::vec4 p2);
	Triangle();
};

class Plane
{
public:
	glm::vec4 normal;
	glm::vec4 point;
	unsigned int uSurf;
	unsigned int reserved0;
	unsigned int reserved1;
	unsigned int reserved2;

	Plane(glm::vec4 normal, glm::vec4 point);
	Plane();
};

class GeometryDescriptor{
public:
	int numSpheres;
	int numPlanes;
	int numTriangles;

	GeometryDescriptor( int a, int b, int c );
};

class Surface
{
public:
	glm::vec4 vColor;
	glm::vec4 vEmissive;   
	Surface(glm::vec4 color, glm::vec4 emissive);
	Surface();
};

#endif
