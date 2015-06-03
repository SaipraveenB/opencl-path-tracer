#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

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
  int iSphere;

} Intersection;

// Returns a normalized vector corresponding to the camera context and screen position.
float4 shoot_ray( float2 screen, __global Camera* cam )
{

  float z = 1.0f; // Assume unit Z Depth.

  float4 vCompY = cam->vLookAt * z + cam->vUp * z * screen.y;
  float4 vCompX = cam->vLookAt * z + cross( cam->vUp, cam->vLookAt ) * z * screen.x;

  //return (screen.x, screen.y, 1.0f, 0.0f );
  return normalize( vCompX + vCompY );

}


float getDiscriminant(float4 lStart, float4 lDir, Sphere sphere)
{
    float dot1 = dot(lDir, (lDir - sphere.vPos));
    dot1 *= dot1;
    float dist = dot(lDir - sphere.vPos,lDir - sphere.vPos);
    return dot1 - dist + (sphere.fRadius*sphere.fRadius);
}

Intersection ray_intersect( float4 vPos, float4 vDir, __global Sphere* pSphere, int maxSpheres )
{
    Intersection vIntPoint;
    vIntPoint.iSphere = -1;
    vIntPoint.vPos = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

    // AAAARGGGHHH... Y U name this INT instead of FLOAT? Took 3 HOURS to figure it out. :P
    float curDisc = 0.0f;
    //minDist = (float)1000000.0f;
    float minDist = 10000.0f;

    /*for( int i = 0; i < maxSpheres; i++ ){
      curDisc ++;
    }*/

    //float distA, distB;
    for(int i=0;i < maxSpheres;i++)
    {
        //curDisc = getDiscriminant(vPos,vDir,pSphere[i]);

        float dot1 = dot(vDir, (vPos - pSphere[i].vPos));
        dot1 = dot1 * dot1;
        float dist = dot(vPos - pSphere[i].vPos,vPos - pSphere[i].vPos);
        curDisc = dot1 - dist + (pSphere[i].fRadius*pSphere[i].fRadius);

        //curDisc = -1.0f;
        if(curDisc < 0)
            continue;

        //i = 1;
        float distA = -(dot(vDir,vPos - pSphere[i].vPos)) - curDisc;
        float distB = distA + 2*curDisc;
        //float distA = 1.0f;
        //float distB = 1.5f;
        if(distA < 0)
            distA = minDist;
        if(distB < 0)
            distB = minDist;

        //curDisc = 1;
        //minDist = 1.0f;
        //vIntPoint.iSphere = 0;


        if(distA < minDist)
        {
            vIntPoint.iSphere = i;
            vIntPoint.vPos = vPos + distA * vDir;
            minDist = distA;
        }

        if(distB < minDist)
        {
            minDist = distB;
            vIntPoint.iSphere = i;
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

   float xt = (xf*2.0 - 1.0f) * (1920.0f / 1080.0f);
   float yt = yf*2.0 - 1.0f;

   float2 screen = (float2)( xt, yt );

   // Shoot ray through screen.
   float4 ray = (shoot_ray( screen, pCamera ));
   //float4 ray = cross( pCamera[0].vUp, pCamera[0].vLookAt );
   //ray = ray * sign( ray );
   //ray = pCamera[0].vUp;

   Intersection patch = ray_intersect( pCamera[0].vPos, ray, pObjects, 1 );

   float mag = 0.0f;
   if( patch.iSphere != -1 ){
    mag = 1.0f;
   }
   // Note: no samplers have been used so indexing is directly through integers.
   write_imagef( imgOut, (int2)(x,y), (float4)( mag, mag, mag, 1.0f ) );

}
