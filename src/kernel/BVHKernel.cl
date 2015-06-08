//Morton Code calculation
//Fast Parallel Radix sort on GPU.
//http://amd-dev.wpengine.netdna-cdn.com/wordpress/media/2013/01/AMD_OpenCL_Tutorial_SAAHPC2010.pdf

//Assumed a primitive called Object that can be used to represent any object.

typedef struct __BVol Volume;
struct Volume
{
  Volume* left;
  Volume* right;
  Object myObject;
  float2 xBounds;
  float2 yBounds;
  float2 zBounds;
};

typedef struct __Int
{
  Object primObject;
  int primType;
  float4 intPoint;
} Intersection;     

//If the given ray intersects the given bounding volume cube.
bool doesIntersect(Volume* vol, float4 vPos, float4 vDir)
{
  
}

//Returns the Intersection pt of ray and the first closest  bounding volume.
Intersection getIntersection(Volume* vol, float4 vPos, float4 vDir)
{
  if(vol == NULL)
  {
    Intersection noInt;
    noInt.primType = -1;
    return noInt;
  }
}
