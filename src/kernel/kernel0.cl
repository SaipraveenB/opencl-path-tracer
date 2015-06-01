#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

__constant sampler_t sampler =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

typedef struct _Camera
{
  float3 vPos;
  float3 vLookAt;
  float3 vUp;

} Camera;

// Primitive types and listings
typedef struct _Sphere
{

  float3 vPos;
  float fRadius;

} Sphere;

typedef struct _Intersection{

  float3 vPos;
  Sphere* pSphere;

} Intersection;

// Returns a normalized vector corresponding to the camera context and screen position.
float3 shoot_ray( float2 screen, Camera* cam )
{

  float z = 1.0f; // Assume unit Z Depth.

  float3 vCompY = vLookAt * z + vUp * z * screen.y;
  float3 vCompX = vLookAt * z + cross( vUp, vLookAt ) * screen.x;

  return normalize( vCompX + vCompY );
}


inline float getDiscriminant(float3 lStart, float3 lDir, Sphere sphere)
{
    float dot1 = dot(lDir, (lDir - sphere.vPos));
    dot1 *= dot1;
    float dist = dot(lDir - sphere.vPos,lDir - sphere.vPos);
    return dot1 - dist + (sphere.fRadius*sphere.fRadius);
}

Intersection ray_intersect( float3 vPos, float3 vDir, Sphere* pSphere, int maxSpheres )
{
    Intersection vIntPoint;
    vIntPoint.pSphere = NULL;
    int curDisc;
    float minDist = 340282346638528859811704183484516925440.0;
    float distA, distB;
    for(int i=0;i < maxSpheres;i++)
    {
        curDisc = getDiscriminant(vPos,vDir,pShere[i]);
        if(curDisc < 0)
            continue;
        distA = -(dot(vDir,vPos - pShere.vPos)) - curDisc;
        distB = distA + 2*curDisc;
        if(distA < 0)
            distA = minDist;
        if(distB < 0)
            distB = minDist;
        if(distA < minDist)
        {
            minDist = distA;
            vIntPoint.pSphere = pSphere + i;
            vIntPoint.vPos = vPos + distA * vDir;
        }
        if(distB < minDist)
        {
            minDist = distB;
            vIntPoint.pSphere = pSphere + i;
            vIntPoint.vPos = vPos + distB * vDir;
        }
    }

    return vIntPoint;
}

__kernel void path_trace( __global Sphere* pObjects, __global Camera* pCamera, __write_only image2d_t imgOut )
{
   size_t x = get_global_id(0);
   size_t y = get_global_id(1);

  // Normalized coords.
   float xf = ((float)x)/1920.0f;
   float yf = ((float)y)/1080.0f;
   float2 screen = ( xf, yf );

   // Shoot ray through screen.
   float3 ray = shoot_ray( screen, pCamera );

   float xt = ( xf*2.0 - 1.0f );
   float yt = yf*2.0 - 1.0f;

   // Note: no samplers have been used so indexing is directly through integers.
   write_imagef( imgOut, (int2)(x,y), (float4)( ray.x , ray.y, 0.0f, 1.0f ) );

}
