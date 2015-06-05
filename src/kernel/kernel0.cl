#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define PIX_X 0
#define PIX_Y 0

#define log(x) if( PIX_X == get_global_id(0) && PIX_Y == get_global_id(1) ) printf(x);

__constant sampler_t sampler =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

typedef struct _Camera
{

  float4 vPos;
  float4 vLookAt;
  float4 vUp;

} Camera;

// Primitive types and listings
typedef struct _Sphere
{

  float4 vPos;
  float fRadius;

} Sphere;

typedef struct _Intersection{

  float4 vPos;
  float4 vDir;
  uint bInter;

} Intersection;

typedef struct _FullIntersection{
  Intersection kInter;
  uint uPrimType;
  int uIndex;
  float4 vDir;
  // Put other stuff like BRDF, normals, texture coords etc. here,
  // to be calculated by ray_intersect after determining the hit.

} FullIntersection;

typedef struct _GeometryDescriptor{
  int numSpheres;
  int numTriangles;
  int numPlanes;
} GeometryDescriptor;

typedef struct __Triangle
{
    //CCW naming. Right bottom vertex is a.
    float4 a;
    float4 b;
    float4 c;
} Triangle;

typedef struct __Plane
{
    float4 normal;
    float4 pt;
} Plane;



// Returns a normalized vector corresponding to the camera context and screen position.
float4 shoot_ray( float2 screen, __global Camera* cam )
{

  float z = 1.0f; // Assume unit Z Depth.

  float4 vCompY = cam->vLookAt * z + cam->vUp * z * screen.y;
  float4 vCompX = cam->vLookAt * z + cross( cam->vUp, cam->vLookAt ) * z * screen.x;

  //return (screen.x, screen.y, 1.0f, 0.0f );
  return normalize( vCompX + vCompY );

}

float get_discriminant( float4 lStart, float4 lDir, __global Sphere* sphere, int i)
{
    float dot1 = dot(lDir, (lStart - sphere[i].vPos));
    dot1 = dot1 * dot1;
    float dist = dot(lStart - sphere[i].vPos,lStart - sphere[i].vPos);
    return dot1 - dist + (sphere[i].fRadius*sphere[i].fRadius);
}

float4 getNormal( __global Sphere* sphere, int i, float4 vPos )
{
    return normalize( vPos - sphere[i].vPos);
}

Intersection rayPlaneIntersect( float4 vPos, float4 vDir, __global Plane* plane, int i )
{
    Intersection x;
    x.vPos = (float4)(0,0,0,0);
    x.vDir = (float4)(0,0,0,0);

    int angle = dot(vDir, plane[i].normal);
    if(angle == 0)
        return x;

    x.vPos = dot(plane[i].pt - vPos, plane[i].normal);
    x.vDir = plane[i].normal;

    return x;
}

Intersection rayTriangleIntersect(  float4 vPos, float4 vDir, __global Triangle* delta, int i )
{
    Intersection intPoint;
    intPoint.bInter = 0;
    float4 normal = normalize(cross(delta[i].b - delta[i].a, delta[i].c - delta[i].a));
    float val = dot(vDir,normal);

    if(val == 0)
        return intPoint;

    float lineDist = dot(delta[i].a - vPos,normal)/val;
    float4 point = vPos + lineDist*vDir;

    //Check if the point lies inside the triangle.
    float dotA = dot(cross(point-delta[i].a,delta[i].c-delta[i].a),normal);
    float dotB = dot(cross(point-delta[i].c,delta[i].b-delta[i].c),normal);
    float dotC = dot(cross(point-delta[i].b,delta[i].a-delta[i].b),normal);

    //Either of the 3 dot products is negative if point is outside traingle.
    if(dotA <0 || dotB < 0 || dotC < 0)
        return intPoint;

    intPoint.bInter = i;
    intPoint.vPos = point;

    return intPoint;
}

Intersection raySphereIntersect( float4 vPos, float4 vDir, __global Sphere* pSphere, int i ){

  Intersection vIntPoint;
  vIntPoint.bInter = 0;
  vIntPoint.vPos = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

  float curDisc = get_discriminant(vPos,vDir,pSphere,i);

  if(curDisc < 0)
      return vIntPoint;

  //i = 1;
  float distA = -(dot(vDir,vPos - pSphere[i].vPos)) - curDisc;
  float distB = distA + 2*curDisc;

  if( distA < 0 && distB < 0 )
      return vIntPoint;


  if( distB < 0 || distA < distB)
  {
      //minDist = distA;
      vIntPoint.bInter = 1;
      vIntPoint.vPos = vPos + distA * vDir;
      //vIntPoint.vDir = getNormal( pSphere, i, vIntPoint.vPos );
  }

  if( distA < 0 || distB < distA)
  {
      //minDist = distB;
      vIntPoint.bInter = 1;
      vIntPoint.vPos = vPos + distB * vDir;
      //vIntPoint.vDir = getNormal( pSphere, i, vIntPoint.vPos );
  }

  return vIntPoint;
}

// Intersect with all types of objects.
// Primitive form of ray intersection: brute force.
FullIntersection ray_intersect( float4 vPos, float4 vDir, __global Sphere* pSpheres, __global Plane* pPlanes, __global Triangle* pTriangles, __global GeometryDescriptor* pDescriptor )
{

    //Intersection vIntPoint;
    //vIntPoint.bInter = -1;
    //vIntPoint.vPos = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

    float curDisc = 0.0f;
    //minDist = (float)1000000.0f;
    float minDist = 10000.0f;

    Intersection inter;

    //float minDist = 10000.0f;

    float primType = -1;
    float index = -1;

    //float distA, distB;
    for(int i=0;i < pDescriptor[0].numSpheres;i++)
    {
      inter = raySphereIntersect( vPos, vDir, pSpheres, i );
      //printf("%d, %d: Sphere: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
      if( inter.bInter != 0 ){
        float dist = dot( inter.vPos, vDir );
        if( dist < minDist ){
          minDist = dist;
          primType = 0;
          index = i;
        }
      }

    }

    for(int i=0;i < pDescriptor[0].numPlanes;i++)
    {
      inter = rayPlaneIntersect( vPos, vDir, pPlanes, i );
      //printf("%d, %d: Plane: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
      if( inter.bInter != 0 ){
        float dist = dot( inter.vPos, vDir );
        if( dist < minDist ){
          minDist = dist;
          primType = 1;
          index = i;
        }
      }
    }

    for(int i=0;i < pDescriptor[0].numTriangles;i++)
    {
      inter = rayTriangleIntersect( vPos, vDir, pTriangles, i );
      //printf("%d, %d: Triangle: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
      if( inter.bInter != 0 ){
        float dist = dot( inter.vPos, vDir );
        if( dist < minDist ){
          minDist = dist;
          primType = 2;
          index = i;
        }
      }
    }

    FullIntersection full;

    full.uPrimType = primType;
    full.uIndex = index;
    full.kInter = inter;

    // To find the normal.
    if( full.uPrimType == 0 ){
      full.vDir = normalize( pSpheres[full.uIndex].vPos - full.kInter.vPos );
    }else if( full.uPrimType == 1){
      full.vDir = pPlanes[full.uIndex].normal;
    }else if( full.uPrimType == 2){
      int index = full.uIndex;
      full.vDir = normalize( cross( pTriangles[index].a - pTriangles[index].b, pTriangles[index].b - pTriangles[index].c ) );
    }

    //Flip normal if necessary.
    if( dot( vDir, full.vDir ) > 0 ){
      full.vDir = -full.vDir;
    }

    return full;
}

// TODO: Need another kernel for quickly building octree out of Sphere, Triangle and Plane.
/*
  As of now the path tracer supports 3 primitives.
  Spheres, Planes and Triangles.
  Triangles are extremely versatile as they can be made into
*/
__kernel void path_trace( __global Camera* pCamera,
                          __write_only image2d_t imgOut,
                          __read_only image2d_t imgIn,
                          __global Sphere* pSpheres,
                          __global Triangle* pTriangles,
                          __global Plane* pPlanes,
                          __global GeometryDescriptor* pDesc )
{
   size_t x = get_global_id(0);
   size_t y = get_global_id(1);
   if( x == 0 && y == 0 )
    printf("X,Y: %d, %d\n", x, y);
  // Normalized coords.
   float xf = ((float)x)/320.0f;
   float yf = ((float)y)/240.0f;

   float xt = (xf*2.0 - 1.0f) * (320.0f / 240.0f);
   float yt = yf*2.0 - 1.0f;

   float2 screen = (float2)( xt, yt );

   float4 ray = (shoot_ray( screen, pCamera ));

   FullIntersection patch = ray_intersect( pCamera[0].vPos, ray, pSpheres, pPlanes, pTriangles, pDesc );

   float4 col;
   if( patch.kInter.bInter ){
     //float mag = distance( patch.kInter.vPos, pCamera[0].vPos )/10.0f;
     //col = (float4)(mag, mag, mag, 1.0f);
    if( patch.uPrimType == 0 )
        col = (float4)( 0.0f, 1.0f, 0.0f, 0.0f );
    else if( patch.uPrimType == 1 )
        col = (float4)( 1.0f, 0.0f, 0.0f, 0.0f );
    else if( patch.uPrimType == 2 )
        col = (float4)( 0.0f, 0.0f, 1.0f, 0.0f );

   }else{
     col = (float4)( 1,0,0,1 );
   }

   col = (float4)( pDesc[0].numSpheres, pDesc[0].numPlanes, pDesc[0].numTriangles, 1.0f );
   float mag = 0.0f;

   /*if( patch.bInter != -1 ){
     float4 vec = patch.vPos - pObjects[patch.bInter].vPos;

     float mag = -dot( vec, pCamera[0].vLookAt );
     //col = (float4)( vec.x * vec.x, vec.y * vec.y, vec.z * vec.z, 1.0f );
     col = (float4)( mag, mag, mag, 1.0f );
   }else{
     col = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );
   }*/
   // Note: no samplers have been used so indexing is directly through integers.
   write_imagef( imgOut, (int2)(x,y), col );

}
