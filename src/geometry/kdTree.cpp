#include "kdTree.h"
#include <stdgl.h>
#include <primitives.h>

#define SQR (a) a*a

int BBox::longestAxis()
{
        float xLength = xBounds.y - xBounds.x;
        float yLength = yBounds.y - yBounds.x;
        float zLength = zBounds.y - zBounds.x;

        //Find max;
        float m = a;
        (m < b) && (m = b); //these are not conditional statements.
        (m < c) && (m = c); //these are just boolean expressions.

        return (m == a ? 1 : m == b ? 2 : 3);
}

//Check if tri/sphere intersects the box / is member of box.
bool isMember(BBox &box, Triangle* tri)
{
    //Caching to make access faster.
    glm::vec4 p0(tri->p0);
    glm::vec4 p1(tri->p1);
    glm::vec4 p2(tri->p2);
    bool result0 = (p0.x > box.xBounds.x && p0.x < box.xBounds.y) && (p0.y > box.yBounds.x && p0.y < box.yBounds.y) && (p0.z > box.zBounds.x && p0.z < box.zBounds.y);
    bool result1 = (p1.x > box.xBounds.x && p1.x < box.xBounds.y) && (p1.y > box.yBounds.x && p1.y < box.yBounds.y) && (p1.z > box.zBounds.x && p1.z < box.zBounds.y);
    bool result2 = (p2.x > box.xBounds.x && p2.x < box.xBounds.y) && (p2.y > box.yBounds.x && p2.y < box.yBounds.y) && (p2.z > box.zBounds.x && p2.z < box.zBounds.y);
    return (result0 || result1 || result2);
}

bool isMember(BBox &box, Sphere* sphere)
{
    //Caching to make access faster.
    glm::vec4 center = sphere->center;
    float rad_sqr = sphere->radius * sphere->radius;
    glm::vec3 minCorner(box.xBounds.x,box.yBounds.x,box.zBounds.x);
    glm::vec3 maxCorner(box.xBounds.y,box.yBounds.y,box.zBounds.y);
    float dmin = 0.0f;

    ( center.x < minCorner.x ) && (dmin += SQR( center.x - minCorner.x ));
    ( center.x > maxCorner.x ) && (dmin += SQR( center.x - maxCorner.x ));
    ( center.y < minCorner.y ) && (dmin += SQR( center.y - minCorner.y ));
    ( center.y > maxCorner.y ) && (dmin += SQR( center.y - maxCorner.y ));
    ( center.z < minCorner.z ) && (dmin += SQR( center.z - minCorner.z ));
    ( center.z > maxCorner.z ) && (dmin += SQR( center.z - maxCorner.z ));

    return dmin <= rad_sqr;
}

//Splits the BBox into two BBox objects using "mid" as value for splitting in "axis" axis.
void BBox::split(float mid,int axis, BBox &left, BBox &right)
{
    left.xBounds = this->xBounds;
    left.yBounds = this->yBounds;
    left.zBounds = this->zBounds;
    right.xBounds = this->xBounds;
    right.yBounds = this->yBounds;
    right.zBounds = this->zBounds;

    (axis == 1) && (left.xBounds.y = mid, right.xBounds.x);
    (axis == 2) && (left.yBounds.y = mid, right.yBounds.x);
    (axis == 3) && (left.zBounds.y = mid, right.yBounds.x);

    int size = this->tris.size();
    for(int index = 0;index < size;index++)
        isMember(left,left.tris[i]) ? left.tris.push_back(this->tris[i]) : right.tris.push_back(this->tris[i]);

    size = this->sphres.size();
    for(int index = 0;index < size;index++)
        isMember(left,left.tris[i]) ? left.tris.push_back(this->tris[i]) : right.tris.push_back(this->tris[i]);

    return;
}

//Required to find intersection point.
bool BBox::doesIntersect(glm::vec3 vPos, glm::vec3 vDir)
{
    float tNear = INT_MIN,tFar = INT_MAX;
    float curNear,curFar,tempNear,tempFar;
    bool result = true;

    //TODO :: Do the following for each pair of opp faces.
    //------------------------------
    glm::vec3 p1(this->xBounds.x,this->yBounds.y, this->zBounds.x);
    glm::vec3 p2(this->xBounds.x,this->yBounds.x, this->zBounds.x);
    glm::vec3 p3(this->xBounds.y,this->yBounds.x, this->zBounds.x);
    //On the opposite face. Only one point required - normal is the same.
    glm::vec3 p4(this->xBounds.x,this->yBounds.y, this->zBounds.x);

    glm::vec3 norm = glm::fastNormalize(glm::cross(p1-p2,p3-p2));

    float sec = 1/glm::dot(norm,vDir);
    (sec >= 0) && (result = false);
    tempNear = glm::dot(p2-vPos,norm)*sec;
    tempFar = glm::dot(p4-vPos,norm)*sec;
    curNear = tempNear;
    curFar = tempFar;

    (tempNear > tempFar) && (curNear = tempFar,curFar = tempNear);
    (curNear > tNear) && (tNear = curNear);
    (curFar < tFar) && (tFar = curFar);
    //--------------------------------------

    return (tNear  > tFar);
}

//Required to stop the k-D Tree traversal.
bool BBox::stop()
{
    return (this->tris.size() + this->spheres.size() <= KD_TREE_MAX_THRESHOLD);
}

//TODO :: Check whether this is needed.
//Returns the intersection point if it exists. Returns (-2,-2,-2) otherwise.
glm::vec3 getX(Triangle *tri, glm::vec2 vPos, glm::vec2 vDir)
{

}
glm::vec3 getX(Sphere *tri, glm::vec2 vPos, glm::vec2 vDir)
{

}
