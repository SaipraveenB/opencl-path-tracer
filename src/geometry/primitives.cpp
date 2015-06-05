#include <geometry/primitives.h>

Sphere::Sphere(glm::vec4 center, float rad)
{

    this->center = center;
    this->radius = rad;

}

Sphere::Sphere(){

    this->center = glm::vec4(0.0f,0.0f,0.0f,0.0f);
    this->radius = 1.0f;

}


Triangle::Triangle(glm::vec4 p0, glm::vec4 p1, glm::vec4 p2 )
{

    this->p0 = p0;
    this->p1 = p1;
    this->p2 = p2;

}

Triangle::Triangle()
{

    this->p0 = glm::vec4(0.0f,0.0f,0.0f,0.0f);;
    this->p1 = glm::vec4(0.0f,0.0f,0.0f,0.0f);;
    this->p2 = glm::vec4(0.0f,0.0f,0.0f,0.0f);;

}

Plane::Plane(glm::vec4 normal, glm::vec4 point){

    this->normal = normal;
    this->point = point;

}

Plane::Plane(){

    this->normal = glm::vec4(0.0f,0.0f,0.0f,0.0f);
    this->point = glm::vec4(0.0f,0.0f,0.0f,0.0f);

}


GeometryDescriptor::GeometryDescriptor( int a, int b, int c ){

    this->numSpheres = a;
    this->numTriangles = b;
    this->numPlanes = c;

}
