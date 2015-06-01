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

Intersection ray_intersect( float3 vPos, float3 vDir, Sphere* pSphere )
{
  // For now calculate line distance to each sphere and get the nearest intersection.
}


__kernel void path_trace( __global Sphere* pObjects, __global Camera* pCamera, __write_only image2d_t imgOut )
{
   size_t x = get_global_id(0);
   size_t y = get_global_id(1);

  // Normalized coords.
   float xf = ((float)x)/1920.0f;
   float yf = ((float)y)/1080.0f;
   float2 screen = ( xf, yf );

   shoot_ray( screen, pCamera );

   float xt = ( xf*2.0 - 1.0f );
   float yt = yf*2.0 - 1.0f;

   // Note: no samplers have been used so indexing is directly through integers.
   write_imagef( imgOut, (int2)(x,y), (float4)( xf , yf, 0.0f, 1.0f ) );

}
