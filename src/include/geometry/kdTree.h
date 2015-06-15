#ifndef KDTREE_H
#define KDTREE_H

#include <stdgl.h>
#include <primitives.h>

#define KD_TREE_MAX_THRESHOLD 5

class BBox
{
public:
    glm::vec2           xBounds;
    glm::vec2           yBounds;
    glm::vec2           zBounds;
    vector<Triangle*>   tris;
    vector<Sphere*>     spheres;

    //Returns the axis to split the BBox.
    int longestAxis();

    //Splits the BBox into two BBox's "left" and "right" using "mid" as value for splitting in "axis" axis.
    void split(float mid,int axis, BBox &left, BBox &right);

    //Required to find intersection point.
    bool doesIntersect(glm::vec3 vPos, glm::vec3 vDir);

    //Required to stop the k-D Tree traversal.
    bool stop();
}

//Returns the intersection point if it exists. Returns (-2,-2,-2) otherwise.
glm::vec3 getX(Triangle *tri, glm::vec2 vPos, glm::vec2 vDir);
glm::vec3 getX(Sphere *tri, glm::vec2 vPos, glm::vec2 vDir)
