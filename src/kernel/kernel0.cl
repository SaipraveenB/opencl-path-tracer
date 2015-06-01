#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

__constant sampler_t sampler =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;


typedef struct _Camera{
  float3 vPos;
  float3 vLookAt;
  float3 vUp;
} Camera;

typedef struct _Object{
  uint uType;
  uint uPrimitive;
} Object;

// Primitive types and listings
typedef struct _Sphere{
  float3 vPos;
} Sphere;

// Returns a normalized vector corresponding to the camera context and screen position.
float3 shootVector( int2 screen, Camera* cam ){
  return (float3)(0.0f,0.0f,0.0f);
}

__kernel void square(__global float * out, __global float * in, __write_only image2d_t imgOut) {
   size_t x = get_global_id(0);
   size_t y = get_global_id(1);

  // Normalized coords.
   float xf = ((float)x)/1920.0f;
   float yf = ((float)y)/1080.0f;

   float xt = ( xf*2.0 - 1.0f );
   float yt = yf*2.0 - 1.0f;

   float mag = 1.0f;


   // Note: no samplers have been used so indexing is directly through integers.
   write_imagef( imgOut, (int2)(x,y), (float4)( xf , yf, 0.0f, 1.0f ) );

}
