#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define PIX_X 533
#define PIX_Y 460

#define WIDTH 1067
#define HEIGHT 600
#define SAMPLES 10


//#define CMD_DEBUG
//#define INDIRECT_ONLY

#define BLUR_FACTOR 2.0f
#ifdef CMD_DEBUG
#define cmdlog(x, ...) if( PIX_X == get_global_id(0) && PIX_Y == get_global_id(1) ) printf(x, __VA_ARGS__);
#else
#define cmdlog(x, ...) ;
#endif

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
  uint uSurface;
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
  float4 vDiffuse;
  float4 vEmissive;

} FullIntersection;

typedef struct _GeometryDescriptor{
  int numSpheres;
  int numPlanes;
  int numTriangles;
} GeometryDescriptor;

typedef struct _ImageDescriptor{
  int numSamples;
  int bChanged;
  int sampleRate
} ImageDescriptor;

typedef struct __Triangle
{
    //CCW naming. Right bottom vertex is a.
    float4 a;
    float4 b;
    float4 c;
    uint uSurface;
} Triangle;

typedef struct __Plane
{
    float4 normal;
    float4 pt;
    uint uSurface;
} Plane;


// Surface properties container.
typedef struct _Surface
{
    float4 vColor;
    float4 vEmissive;
} Surface;

// CORE RANDOMIZERS

// XORSHIFT128 RNG for decent randomness.
/*uint rng_next( __global uint4* state, int k ) {
    cmdlog("State: %d,%d,%d,%d\n", (int)state[k].x, (int)state[k].y, (int)state[k].z,(int)state[k].w );
    uint t = state[k].x ^ ( state[k].x << 11);
    state[k].x = state[k].y; state[k].y = state[k].z; state[k].z = state[k].w;
    state[k].w = state[k].w ^ (state[k].w >> 19) ^ t ^ (t >> 8);
    return state[k].w;
}*/
uint rng_next( private uint4* state, int k ) {
    state[0].w = 48271*state[0].w;
    return state[0].w;
}
float get_rng_float( private uint4* state, int k ){
    return (float)rng_next( state, k )/(float)UINT_MAX;
}

// UTIL METHODS FOR RANDOM SAMPLING

// Uniformly sample a hemisphere oriented toward a given unit vector.
float4 random_sample_hemisphere( private uint4* state, int k, float4 vDir ){
    //uint4 ray = ( rng_next( state, k ), rng_next( state, k ) ,rng_next( state, k ), rng_next( state, k ) );
    //vDir.w = 0.0f;

    float p0 = get_rng_float( state, k );
    float p1 = get_rng_float( state, k );

    float p2 = get_rng_float( state, k );
    // Form perpendicular vector.

    float4 vecRand = (float4) ( p2*2.0 - 1.0f, 0.0f, 1.0f, 0.0f );
    //float4 vPerp = (float4) ( vDir.y, -vDir.x, 0.0f, 0.0f );
    //if( dot( vecRand, vDir ) == 0 )
        //vecRand = (float4) (1.0f,0.0f,0.0f, 0.0f);

    float4 vPerp = fast_normalize( cross( vDir, vecRand ) );
    // Form the third perpendicular vector for an orthonormal set.
    float4 vPerp2 = cross( vDir, vPerp );

    cmdlog("vPerp: %f, %f, %f, %f\n", vPerp.x, vPerp.y, vPerp.z, vPerp.w);
    cmdlog("vPerp2: %f, %f, %f, %f\n", vPerp2.x, vPerp2.y, vPerp2.z, vPerp2.w);
    // treat p0 as lift angle and p1 as sweep angle.

    // trasform p1 from (0,1) to (-1,+1)
    p1 = p1*6.283185307f;
    p0 = p0*1.570796326f;
    cmdlog("p0,p1: %f,%f\n", p0, p1);
    /*float cos0;
    float sin0 = sincos( p0, &cos0 );
    
    float cos1;
    float sin1 = sincos( p1, &cos1 );*/
    
    float4 ray = ( native_cos( p0 ) * vDir ) + ( native_sin( p0 ) * ( native_cos( p1 ) * vPerp + native_sin( p1 ) * vPerp2 ) );
    //ray.w = 0.0f;
    return ray;
}

float4 random_sample_cone( private uint4* state, int k, float4 vDir, float halfAngle ){
    
    float p0 = get_rng_float( state, k );
    float p1 = get_rng_float( state, k );

    float p2 = get_rng_float( state, k );
    // Form perpendicular vector.

    float4 vecRand = (float4) ( p2*2.0 - 1.0f, 0.0f, 1.0f, 0.0f );
    //float4 vPerp = (float4) ( vDir.y, -vDir.x, 0.0f, 0.0f );
    //if( dot( vecRand, vDir ) == 0 )
        //vecRand = (float4) (1.0f,0.0f,0.0f, 0.0f);

    float4 vPerp = fast_normalize( cross( vDir, vecRand ) );
    // Form the third perpendicular vector for an orthonormal set. 	
    float4 vPerp2 = cross( vDir, vPerp );

    //cmdlog("vPerp: %f, %f, %f, %f\n", vPerp.x, vPerp.y, vPerp.z, vPerp.w);
    //cmdlog("vPerp2: %f, %f, %f, %f\n", vPerp2.x, vPerp2.y, vPerp2.z, vPerp2.w);
    // treat p0 as lift angle and p1 as sweep angle.

    // trasform p1 from (0,1) to (-1,+1)
    p1 = p1*6.283185307f;
    p0 = ( 1.0f - ( acos(p0) * 2.0f * M_1_PI ) ) * halfAngle;
    //p0 = half_sqrt( p0 ) * halfAngle;
    //cmdlog("p0,p1: %f,%f\n", p0, p1);

    float4 ray = ( native_cos( p0 ) * vDir ) + ( native_sin( p0 ) * ( native_cos( p1 ) * vPerp + native_sin( p1 ) * vPerp2 ) );
    //ray.w = 0.0f;
    return ray;
}



//Returns a float2 coordinate to shoot a ray.
float2 findRandomPoint(private uint4* state, int k, float2 pix)
{
    cmdlog("Random point:  %f, %f", pix.x, pix.y);
    float2 newPt;
    float scaleX = get_rng_float(state,k);
    float scaleY = get_rng_float(state,k);
    float deltaX = native_divide( BLUR_FACTOR, (float)WIDTH );
    float deltaY =  native_divide( BLUR_FACTOR, (float)HEIGHT );
    newPt.x = pix.x + scaleX*deltaX;
    newPt.y = pix.y + scaleY*deltaY;
    cmdlog(":%f, %f\n", newPt.x, newPt.y);
    return newPt;
}


// Returns a normalized vector corresponding to the camera context and screen position.
float4 shoot_ray( float2 screen, __global Camera* cam )
{

  float z = 1.0f; // Assume unit Z Depth.

  float4 vCompY = cam->vLookAt * z + cam->vUp * z * screen.y;
  float4 vCompX = cam->vLookAt * z + cross( cam->vUp, cam->vLookAt ) * z * screen.x;

  //return (screen.x, screen.y, 1.0f, 0.0f );
  return fast_normalize( vCompX + vCompY );

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
    return fast_normalize( vPos - sphere[i].vPos);
}

Intersection rayPlaneIntersect( float4 vPos, float4 vDir, __global Plane* plane, int i )
{
    /*cmdlog("\nPlaneInt\n", 0);
    cmdlog("vPos %f, %f, %f\n", vPos.x, vPos.y, vPos.z);
    cmdlog("vDir %f, %f, %f\n", vDir.x, vDir.y, vDir.z);
    cmdlog("i %d\n", i);
    cmdlog("plane.normal %f, %f, %f\n", plane[i].normal.x, plane[i].normal.y, plane[i].normal.z);
	cmdlog("plane.pt %f, %f, %f\n", plane[i].pt.x, plane[i].pt.y, plane[i].pt.z);
	cmdlog("plane.uSurf %d\n", (int)plane[i].uSurface);*/

    private Plane pl = plane[i];

    Intersection x;
    x.vPos = (float4)(0,0,0,0);
    x.vDir = (float4)(0,0,0,0);
    x.bInter = 0;

    float angle = dot(vDir, pl.normal);
	
	// Enable backface culling.. necesary for creating an enclosure using planes.
    if(angle >= 0)
        return x;
	

	float temp = dot(pl.pt - vPos, pl.normal);
    float desc = dot(temp*pl.normal, vDir);
    if( desc > 0 ){
		
        x.bInter = 1;
	}
    
    float k = native_divide( temp, angle );
    x.vPos = k * vDir + vPos;
    x.vDir = pl.normal;
	
	
    //cmdlog("intPoint.vPos %f, %f, %f\n", x.vPos.x, x.vPos.y, x.vPos.z);
    //cmdlog("intPoint.bInter %d\n", x.bInter);
    return x;
}

Intersection rayTriangleIntersect(  float4 vPos, float4 vDir, __global Triangle* delta, int i )
{
    /*cmdlog("\nTriangleInt\n", 0);
    cmdlog("vPos %f, %f, %f\n", vPos.x, vPos.y, vPos.z);
    cmdlog("vDir %f, %f, %f\n", vDir.x, vDir.y, vDir.z);
    cmdlog("i %d\n", i);
    cmdlog("delta.a %f, %f, %f\n", delta[i].a.x, delta[i].a.y, delta[i].a.z);
	cmdlog("delta.b %f, %f, %f\n", delta[i].b.x, delta[i].b.y, delta[i].b.z);
	cmdlog("delta.uSurf %d\n", (int)delta[i].uSurface);*/
	
    Intersection intPoint;
    intPoint.bInter = 0;
    float4 normal = normalize(cross(delta[i].b - delta[i].a, delta[i].c - delta[i].a));
    float val = dot(vDir,normal);

    if(val > 0)
        return intPoint;

    float lineDist = dot(delta[i].a - vPos,normal)/val;
	if( lineDist < 0 )
		return intPoint;
	
    float4 point = vPos + lineDist*vDir;
    //cmdlog("point %f, %f, %f\n", point.x, point.y, point.z);
    //Check if the point lies inside the triangle.
    float dotA = dot(cross(point-delta[i].a,delta[i].c-delta[i].a),normal);
    float dotB = dot(cross(point-delta[i].c,delta[i].b-delta[i].c),normal);
    float dotC = dot(cross(point-delta[i].b,delta[i].a-delta[i].b),normal);

    //cmdlog("dotA, dotB, dotC %f, %f, %f\n", dotA, dotB, dotC);
    //Either of the 3 dot products is negative if point is outside traingle.
    if(dotA <0 || dotB < 0 || dotC < 0)
        return intPoint;

    intPoint.bInter = 1;
    intPoint.vPos = point;
    //cmdlog("intPoint.vPos %f, %f, %f\n", intPoint.vPos.x, intPoint.vPos.y, intPoint.vPos.z);
    //cmdlog("intPoint.bInter %d\n", intPoint.bInter);
    return intPoint;
}

Intersection raySphereIntersect( float4 vPos, float4 vDir, __global Sphere* pSphere, int i ){

  private Sphere sph = pSphere[i];

  Intersection vIntPoint;
  vIntPoint.bInter = 0;
  vIntPoint.vPos = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

  /*cmdlog("vPos %f, %f, %f\n", vPos.x, vPos.y, vPos.z);
  cmdlog("vDir %f, %f, %f\n", vDir.x, vDir.y, vDir.z);
  cmdlog("i %d\n", i);
  cmdlog("pSphere.vPos %f, %f, %f, %f\n", pSphere[i].vPos.x, pSphere[i].vPos.y, pSphere[i].vPos.z, pSphere[i].vPos.w);
  cmdlog("pSphere.fRadius %f\n", pSphere[i].fRadius);
  cmdlog("pSphere.uSurface %d\n", (int)pSphere[i].uSurface);*/

  //float curDisc = get_discriminant(vPos,vDir,pSphere,i);

  float dot1 = dot(vDir, (vPos - sph.vPos));
  
  if( dot1 > 0 )
    return vIntPoint;

  float dist = dot(vPos - sph.vPos,vPos - sph.vPos);
  float curDisc = dot1 * dot1 - dist + (sph.fRadius * sph.fRadius);
  
  //cmdlog("curDisc: %f\n", curDisc);

  if(curDisc < 0)
      return vIntPoint;

  float discSqrt = half_sqrt(curDisc);

  //i = 1;
  float distA = - dot1 - discSqrt;
  //float distB = distA + 2 * discSqrt;

  if( distA < 0 )
      return vIntPoint;

  vIntPoint.bInter = 1;
  vIntPoint.vPos = vPos + (distA) * vDir;

  //cmdlog("distA: %f, distB: %f\n", distA, distB );
  //cmdlog("intPoint.vPos %f, %f, %f\n", vIntPoint.vPos.x, vIntPoint.vPos.y, vIntPoint.vPos.z);
  //cmdlog("intPoint.bInter %d\n", vIntPoint.bInter);

  return vIntPoint;
}

// Intersect with all types of objects.
// Primitive form of ray intersection: brute force.
FullIntersection ray_intersect( float4 vPos, float4 vDir, __global Sphere* pSpheres, __global Plane* pPlanes, __global Triangle* pTriangles, __global GeometryDescriptor* pDescriptor, __global Surface* pSurfaces , uint originPrim, uint originIndex )
{

    //Intersection vIntPoint;
    //vIntPoint.bInter = -1;
    //vIntPoint.vPos = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

    float curDisc = 0.0f;
    //minDist = (float)1000000.0f;
    float minDist = 10000.0f;

    Intersection inter;

    Intersection currInter;
    Surface surface;
    currInter.bInter = 0;

    //float minDist = 10000.0f;

    float primType = -1;
    float index = -1;

    //float distA, distB;
    for(int i=0;i < pDescriptor[0].numSpheres;i++)
    {
      if( originPrim == 0 && originIndex == (uint)i )
        continue;

      inter = raySphereIntersect( vPos, vDir, pSpheres, i );
      //printf("%d, %d: Sphere: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
      //cmdlog( "sphere: inter with %d: %d\n", i, inter.bInter );
      if( inter.bInter != 0 ){

        float dist = distance( inter.vPos, vPos );
        //cmdlog( "sphere: dist: %f\n", dist  );
        if( dist < minDist ){
          minDist = dist;
          currInter = inter;
          primType = 0;
          index = i;
          surface = pSurfaces[pSpheres[i].uSurface];
        }
      }

    }
	
	//cmdlog("numPLanes: %d\n", pDescriptor[0].numPlanes);
	
    for(int i=0;i < pDescriptor[0].numPlanes;i++)
    {
      if( originPrim == 1 && originIndex == (uint)i )
        continue;
      inter = rayPlaneIntersect( vPos, vDir, pPlanes, i );
      //printf("%d, %d: Plane: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
	  //cmdlog( "plane: inter with %d: %d\n", i, inter.bInter );
      if( inter.bInter != 0 ){

        float dist = distance( inter.vPos, vPos );
		//cmdlog( "plane: dist: %f\n", dist  );
        if( dist < minDist ){
          minDist = dist;
          primType = 1;
          currInter = inter;
          index = i;
          surface = pSurfaces[pPlanes[i].uSurface];
        }
      }
    }

    for(int i=0;i < pDescriptor[0].numTriangles;i++)
    {
      if( originPrim == 2 && originIndex == (uint)i )
        continue;
      inter = rayTriangleIntersect( vPos, vDir, pTriangles, i );
      //printf("%d, %d: Triangle: %d\n", get_global_id(0), get_global_id(1), inter.bInter);
      if( inter.bInter != 0 ){
        float dist = distance( inter.vPos, vPos );
        if( dist < minDist ){
          minDist = dist;
          primType = 2;
          currInter = inter;
          index = i;
          surface = pSurfaces[pTriangles[i].uSurface];
        }
      }
    }

    FullIntersection full;

    full.uPrimType = primType;
    full.uIndex = index;
    full.kInter = currInter;
    full.vDiffuse = surface.vColor;
    full.vEmissive = surface.vEmissive;

    // To find the normal.
    if( full.uPrimType == 0 ){
      full.vDir = fast_normalize( pSpheres[full.uIndex].vPos - full.kInter.vPos );
    }else if( full.uPrimType == 1){
      full.vDir = pPlanes[full.uIndex].normal;
    }else if( full.uPrimType == 2){
      int index = full.uIndex;
      full.vDir = fast_normalize( cross( pTriangles[index].a - pTriangles[index].b, pTriangles[index].b - pTriangles[index].c ) );
    }

    //Flip normal if necessary.
    if( dot( vDir, full.vDir ) > 0 ){
      full.vDir = -full.vDir;
    }

    //cmdlog("FullIntersection: Primtype: %d, index: %d", full.uPrimType, full.uIndex );
    //cmdlog("Intersection: %d\n", full.kInter.bInter);
    return full;
}


// BRDF methods.
// BRDF Methods return a random sampler ray that adheres to the BRDF of the surface.
// Lambertian:
float4 brdf_lambertian( FullIntersection patch, float4 vPos, float4 vDir, private uint4* state, int k ){
    //float p4 = get_rng_float( state, k );
    float4 launchDir = random_sample_cone( state, k, patch.vDir, 3.1415f * 0.50f );
    launchDir.w = 0.5f * M_1_PI;
    //launchDir.w = 1.0f;
    return launchDir;
}

float4 brdf_reflective( FullIntersection patch, float4 vPos, float4 vDir, private uint4* state, int k ){
    float p4 = get_rng_float( state, k );
    
    float4 reflected = ( - 2 * patch.vDir * dot( patch.vDir, vDir ) ) + vDir;
    

    float4 launchDir;
    
    if( p4 > 0.6f )
        launchDir = random_sample_cone( state, k, reflected, 3.1415f * 0.02f );
    else
        launchDir = random_sample_cone( state, k, patch.vDir, 3.1415f * 0.50f );

    launchDir.w = 0.5f * M_1_PI;
    //launchDir.w = 1.0f;
    return launchDir;
}

// Lambertian, Importance sampled.

float4 brdf_lambertian_imp( FullIntersection patch, float4 vPos, float4 vDir, float4 vDirLight, float lightWeight, float halfAngle,private uint4* state, int k ){
    float p4 = get_rng_float( state, k );
   
    float4 dir;
    float angle;
    if( p4 > lightWeight ){
        dir = patch.vDir;
        angle = 0.50f * M_PI;
    }else{
        dir = vDirLight;
        angle = halfAngle;
    }

    float4 launchDir = random_sample_cone( state, k, dir, angle );
    
    float solid = 2 * (float)M_PI * ( 1 - native_cos( halfAngle ) );
    if( solid == 0.0f )
        solid = 0.01f;

    float low = (float)M_1_PI * 0.5f * ( 1 - lightWeight );
    float high = ( 1.0f/solid ) * 0.5f * lightWeight + (float)M_1_PI * 0.5f * ( 1 - lightWeight );
    //float high = (float)M_1_PI * 0.5f * ( 1 - lightWeight );

    //cmdlog("LAMBERT_IMPORTANCE: %f, %f, %f, %f, %f\n", p4, solid, low, high, lightWeight );
    
    if( p4 > lightWeight )
        launchDir.w = ( low ); 
    else
        launchDir.w = ( high );

    return launchDir;
}

// Surface Samplers for lighting.

FullIntersection sampleSphereSurface( __global Surface* pSurface, __global Sphere* pSphere, int i, private uint4* state ){
    FullIntersection fint;

    float p0 = get_rng_float( state, 0 );
    float p1 = get_rng_float( state, 0 );
	
	p0 = p0 * M_PI;
	p0 = acos( 1 - ( 2.0f * M_1_PI * p0 ) );
	
	p1 = p1 * 2.0f * M_PI;

    float4 vDir = (float4)( native_sin( p0 ) * native_sin( p1 ), native_cos( p0 ), native_sin( p0 ) * native_cos( p1 ), 0.0f);
    fint.kInter.vPos = pSphere[i].vPos + pSphere[i].fRadius * vDir;
    fint.kInter.bInter = 1;
    fint.vDir = vDir;

    fint.vEmissive = pSurface[pSphere[i].uSurface].vEmissive;
    
	cmdlog("Light position: %f,%f,%f,%f\n", fint.kInter.vPos.x, fint.kInter.vPos.y, fint.kInter.vPos.z, fint.kInter.vPos.w );
	cmdlog("Light direction: %f,%f,%f,%f\n", fint.vDir.x, fint.vDir.y, fint.vDir.z, fint.vDir.w );
	cmdlog("Light vEmissive: %f,%f,%f,%f\n", fint.vEmissive.x, fint.vEmissive.y, fint.vEmissive.z, fint.vEmissive.w );	
    return fint;
}

__kernel void bi_directional_path_trace( 
                          __global Camera* pCamera,
                          __write_only image2d_t imgOut,
                          __read_only image2d_t imgIn,
                          __global Sphere* pSpheres,
                          __global Plane* pPlanes,
                          __global Triangle* pTriangles,
                          __global GeometryDescriptor* pDesc,
                          __global ImageDescriptor* pImgDesc,
                          __global uint4* randSeed,
                          __global Surface* pSurfaces
                          )
{
    size_t x = get_global_id(0);
    size_t y = get_global_id(1);

    float4 prevValue = read_imagef( imgIn, sampler, (int2)(x,y) );
    int numSamples = pImgDesc[0].numSamples;

    private uint4 state = randSeed[x*HEIGHT + y];

    rng_next( &state, x*HEIGHT + y );
   //if( x == PIX_X && y == PIX_Y )
   /*cmdlog("\nX,Y: %d, %d\n", x, y);
    cmdlog("samples: %d\n", pImgDesc[0].numSamples);
    cmdlog("prevValue: %f, %f, %f, %f\n", prevValue.x, prevValue.y, prevValue.z, prevValue.w);*/
    //cmdlog("sizeof Plane: %d\n", sizeof(Plane) );
	cmdlog("\n\n\n",0);
    // Normalized coords.
    float xf = native_divide( ((float)x), WIDTH*1.0f );
    float yf = native_divide( ((float)y), HEIGHT*1.0f );

    float xt = (xf*2.0 - 1.0f) * ( WIDTH*1.0f / HEIGHT);
    float yt = yf*2.0 - 1.0f;


//    int samples = pImgDesc[0].sampleRate;



    //for( int k = 0; k < samples; k++ ){

    float2 screen = findRandomPoint( &state, x*HEIGHT + y, (float2)( xt, yt ) );

    float4 ray = (shoot_ray( screen, pCamera ));

    FullIntersection patch;// = ray_intersect( pCamera[0].vPos, ray, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, -1, -1 );

    // Second intersection.


    FullIntersection eyePatch;
    eyePatch.vDir = ray;
    eyePatch.kInter.vPos = pCamera[0].vPos;
    eyePatch.uPrimType = 4;
    eyePatch.uIndex = 0;

    float4 col = (float4)(1.0f,1.0f,1.0f,1.0f);
    int i = 0;

    float4 eyeRef = eyePatch.vDir;
    float4 eyePos = eyePatch.kInter.vPos + eyePatch.vDir * 0.00001f;// Push slightly to avoid self-intersection.            
	
    int maxE = 5;
	
	
	FullIntersection currPatchI = eyePatch;
	float4 incomingDirI = glm::vec4( 0, 0, 0, 0 );
	
	float4 eyeVertexList[maxE];
	float4 eyeNormalList[maxE];
		
	cmdlog("\nEye patch trace.\n");
	for( int i = 0; i < maxE; i++ ){
		cmdlog("%d Ray Trace:\n");
		float4 launchRef = currPatchI.vDir;
		float4 launchPos = currPatchI.kInter.vPos + currPatchI.vDir * 0.00001f;// Push slightly to avoid self-intersection.          
		float4 launchDir = brdf_lambertian( currPatchI, (float4)(0,0,0,0), incomingDirI, &state, x*HEIGHT + y );

		
		FullIntersection targetPatch = ray_intersect( launchPos, launchRef, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, currPatchI.uPrimType, currPatchI.uIndex );
		
		if( !targetPatch.kInter.bInter ){
			cmdlog("No intersection\n");
		}
		
		if( targetPatch.vEmissive.w != 0.0f ){
			// Emissive surface. Terminate.
			cmdlog("Target patch is Emissive.")
			break;
		}
		
		eyeVertexList[i] = targetPatch.kInter.vPos;
		eyeNormalList[i] = targetPatch.vDir;		
		
		currPatchI = targetPatch;
		incomingDirI = launchDir;
	}
	
	
   // FullIntersection firstPatch = ray_intersect( eyePos, eyeRef, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, eyePatch.uPrimType, eyePatch.uIndex );

    // Unwrap target point. Soon this needs to be an array.
   // float4 vTestPos = firstPatch.kInter.vPos;
	
    // Handle case where ray strikes emissive surface.
	
	
    if( firstPatch.kInter.bInter )
    {
		cmdlog("First HIT %d, %d\n", firstPatch.uPrimType, firstPatch.uIndex); 
		if( firstPatch.vEmissive.w != 0.0f ){
        	// Break.. we have no further job here.
       		cmdlog("DIRECT HIT LIGHT: %d, %d\n", firstPatch.uPrimType, firstPatch.uIndex); 
        	write_imagef( imgOut, (int2)(x,y), firstPatch.vDiffuse );
		    return;
		}
		
    }else{
		cmdlog("No Intersection: \n", 0);
        write_imagef( imgOut, (int2)(x,y), (float4)(0,0,0,0) );
        return;
    }
	

    float4 incomingDir = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

    
    // Hardcode light values.. soon to be streamed in from host program.
    int lIndex = 0;
    int lType = 0;

    
    // Find the point from light's surface.

    FullIntersection lightPatch = sampleSphereSurface( pSurfaces, pSpheres, lIndex, &state );
    float4 light0Ref = lightPatch.vDir;
    float4 light0Pos = lightPatch.kInter.vPos + lightPatch.vDir * 0.00001f;// Push slightly to avoid self-intersection.            
    //float4 light0Dir;

    // Set initial patch to lightPatch to simulate ray strike from there.
    patch = lightPatch;

    // Assume lambertian BRDF for light patch.
    //float4 light0Dir = brdf_lambertian( patch, (float4)(0,0,0,0), incomingDir, &state, x*HEIGHT + y );

    // Set maximum no of reflections for light ray traversal.
    

    int maxL = 5;


    float pCon = 0.8f;
    // Initialize full contribution to full spectral 0.
    float4 fullContrib = (float4)(0,0,0,0);

    // Accumulator that's diluted every step.
    float4 accContrib = (float4)(1,1,1,1);
	
	float preFactor = (1 - pCon)/(1 - pow( pCon, maxL ));
	int samplesTaken = 0;
	for( int lti = 0; lti < maxL; lti ++ ){
		// Check contribution of current patch to vTestPos.
		
		for( int ei; ei < maxE; ei ++ ){
			cmdlog("\n bounce %d\n", lti);
			
			float4 vTestPos = eyeVertexList[ei];
			float4 vTestDir = eyeNormalList[ei];
						
			float4 testDir = normalize( vTestPos - patch.kInter.vPos );
			cmdlog("testPos: %f, %f, %f, %f\n", vTestPos.x, vTestPos.y, vTestPos.z, vTestPos.w );
			cmdlog("light patch pos: %f, %f, %f, %f\n", patch.kInter.vPos.x, patch.kInter.vPos.y, patch.kInter.vPos.z, patch.kInter.vPos.w );
			cmdlog("ray direction: %f, %f, %f, %f\n", testDir.x, testDir.y, testDir.z, testDir.w );
		
			//cmdlog("test criteria: %d, %d\n", firstPatch.uPrimType, firstPatch.uIndex );
			FullIntersection testPatch = ray_intersect( patch.kInter.vPos, testDir, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, patch.uPrimType, patch.uIndex );
		
			cmdlog("test results: %d, %d\n", testPatch.uPrimType, testPatch.uIndex );
			cmdlog("test results POS: %f, %f, %f, %f\n", testPatch.kInter.vPos.x, testPatch.kInter.vPos.y, testPatch.kInter.vPos.z, testPatch.kInter.vPos.w );
			float4 contrib = (float4)(0,0,0,0);
		
		
			if( firstPatch.uPrimType == testPatch.uPrimType && firstPatch.uIndex == testPatch.uIndex ){
				//if( dot( firstPatch.vDir, testPatch.vDir ) < 0 ){
				if( length( vTestPos - testPatch.kInter.vPos ) < 0.0001f && dot( firstPatch.vDir, testPatch.vDir ) > 0 ){
					cmdlog("LOS acheived.\n", 0);
					float4 ptVec = testPatch.kInter.vPos - patch.kInter.vPos;
					ptVec = ptVec;
				
					float samplingFactor = divide( 1, (float)( ptVec.x*ptVec.x + ptVec.y*ptVec.y + ptVec.z*ptVec.z ) );
					//float samplingFactor = native_divide( 1, 1 + length( ptVec ) );
				
					float r = length( ptVec ) * 1.0f;
					//cmdlog("distance: %f", r);
					//float samplingFactor = 1 - native_divide( r, half_sqrt( r*r + 0.0001f ) );
				
					//if( r < 0.0f )
					//	samplingFactor = 0.0f;
				
				
					//float impFactor = preFactor;
				
					float directionalFactor = dot( vTestDir, testDir ) * dot( patch.vDir, testDir ); 
					directionalFactor *= sign( directionalFactor );
					//if( lti != 0 )
					//	samplingFactor = 1.0f;
					//else 
					//	directionalFactor = 1.0f;
				
					cmdlog("LOS sf: %f, df: %f\n", samplingFactor, directionalFactor);
				
					contrib = accContrib * samplingFactor * directionalFactor; //* native_divide( 1, impFactor );
				
				}
			}
			cmdlog("Adding contrib: %f, %f, %f, %f\n", contrib.x, contrib.y, contrib.z, contrib.w); 
			
			//preFactor *= pCon;
		}
		
		// The 1.0f is the BRDF of the first patch from the eye. we're assuming perfect lambertian. se we set it to 1.0f.
		
		
		// temprarily reject direct light.
#ifdef INDIRECT_ONLY
		if( lti != 0 )
#endif
		fullContrib += contrib;
		
		samplesTaken++;
		
		
		float4 launchRef = patch.vDir;
		float4 launchPos = patch.kInter.vPos + patch.vDir * 0.00001f;// Push slightly to avoid self-intersection.

		float4 launchDir;

		// Soon this part has to be
		
		launchDir = brdf_lambertian( patch, (float4)(0,0,0,0), incomingDir, &state, x*HEIGHT + y );

		//if( lti == 0 )
		//	launchDir = launchRef;
		// Weight at this point.
		float wt = launchDir.w;
		launchDir.w = 0.0f;

		// If the ray is just originating from the sphere.

		if( lti == 0 )
			wt = 0.0f;

		// Contrib so far is unitary.. contribution for the joining ray is computed by using actual BRDF values since it's an operation dependent on 
		// complex geometry and as such can't be efficiently simulated.
		// So random sampling is now uniform and brdf values are used for manipulation.
		

		// roll dice to determine path's fate.
		// But will this bias the estimator?
		float fate = get_rng_float( &state, 0 );
		if( fate > pCon ){
			break;
		}

		// Carry on to another place.
		FullIntersection target = ray_intersect( launchPos, launchDir, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, patch.uPrimType, patch.uIndex );
    
		// Obtain patch information
		if( target.kInter.bInter ){
			if( target.vEmissive.w != 0.0f ){
				//cmdlog( "weight: %f\n", wt );
				//col *= ((float4)patch2.vEmissive/wt) * ( (float)M_1_PI * 0.5f );// sky color.

				cmdlog("%d HIT LIGHT %d,%d\n", i, target.uPrimType, target.uIndex);
				// Break here.. we accidentally seem to have run into a light.
				break;
			}
			//cmdlog( "weight: %f\n", wt );
			//col *= ( (float4)patch2.vDiffuse/wt * (float)M_1_PI * 0.5f );
			cmdlog("%d HIT %d,%d\n", i, target.uPrimType, target.uIndex);
		}else{
			//col *= (float4)(0.3f, 0.3f, 0.4f, 0.3f);// sky color.
			//accContrib *= 0.0f;
			cmdlog("%d Miss\n", i);
			break;
		}
		
		float directionalFactor = -dot( target.vDir, launchDir ) * dot( launchDir, launchRef );
		directionalFactor *= sign( directionalFactor );
		float4 ptVec = ( launchPos - target.kInter.vPos );
		float samplingFactor = native_divide( 1, 1 + (float)( ptVec.x*ptVec.x + ptVec.y*ptVec.y + ptVec.z*ptVec.z ) );
		//float samplingFactor = native_divide( 1, 1 + length( ptVec ) );
		
		//if( lti != 0 ) samplingFactor = 1.0f;
		samplingFactor = 1.0f;
		//else directionalFactor = 1.0f;
		directionalFactor = 1.0f;
		
		cmdlog("sf: %f, df: %f\n", samplingFactor, directionalFactor);		
     	accContrib *= target.vDiffuse * 1.0f * samplingFactor * directionalFactor;
		
		//col = patch.vDiffuse * col;
		
		incomingDir = launchDir;
		patch = target;
    }
    
    
    // Do light_intensity * sigma(contrib) * patch_brdf to get final spectral intensity.
	float4 ptVec = ( firstPatch.kInter.vPos - eyePatch.kInter.vPos );
	float samplingFactor = native_divide( 1, 1 + (float)( ptVec.x*ptVec.x + ptVec.y*ptVec.y + ptVec.z*ptVec.z ) ); 
	samplingFactor = 1.0f;
	if( samplesTaken == 0.0f )
		col = (float)(0,0,0,0);
    else
#ifdef INDIRECT_ONLY
		col = lightPatch.vEmissive * fullContrib * samplingFactor  * native_divide( 1.0f, samplesTaken );// * firstPatch.vDiffuse;
#else
		col = lightPatch.vEmissive * fullContrib * samplingFactor  * native_divide( 1.0f, samplesTaken ) * firstPatch.vDiffuse;
#endif
    //col = (float4) (1.0f, 1.0f, 1.0f, 1.0f);
	
	cmdlog("numSamples: %d\n", numSamples);
    cmdlog("final col: %f, %f, %f, %f; \n", col.x, col.y, col.z, col.w );

    float4 newSample = col;
    float4 currSample = prevValue * native_divide( ((float)numSamples), (numSamples + 1)) + newSample * native_divide( 1.0f, (float)(numSamples + 1) );
    if( numSamples == 0 )
      currSample = newSample;
	
    prevValue = currSample;
    numSamples ++;

    randSeed[x*HEIGHT + y] = state;

    write_imagef( imgOut, (int2)(x,y), currSample );
}

// TODO: Need another kernel for quickly building octree out of Sphere, Triangle and Plane.



/*
  As of now the path tracer supports 3 primitives.
  Spheres, Planes and Triangles.
  Triangles are extremely versatile as they can be made into
*/
/*__kernel void path_trace( __global Camera* pCamera,
                          __write_only image2d_t imgOut,
                          __read_only image2d_t imgIn,
                          __global Sphere* pSpheres,
                          __global Triangle* pTriangles,
                          __global Plane* pPlanes,
                          __global GeometryDescriptor* pDesc,
                          __global ImageDescriptor* pImgDesc,
                          __global uint4* randSeed,
                          __global Surface* pSurfaces
                          )
{
    size_t x = get_global_id(0);
    size_t y = get_global_id(1);

    float4 prevValue = read_imagef( imgIn, sampler, (int2)(x,y) );
    int numSamples = pImgDesc[0].numSamples;

    private uint4 state = randSeed[x*HEIGHT + y];

    rng_next( &state, x*HEIGHT + y );
   //if( x == PIX_X && y == PIX_Y )
    cmdlog("\nX,Y: %d, %d\n", x, y);
    cmdlog("samples: %d\n", pImgDesc[0].numSamples);
    cmdlog("prevValue: %f, %f, %f, %f\n", prevValue.x, prevValue.y, prevValue.z, prevValue.w);
    cmdlog("sizeof Sphere: %d\n", sizeof(Sphere) );
    // Normalized coords.
    float xf = native_divide( ((float)x), WIDTH*1.0f );
    float yf = native_divide( ((float)y), HEIGHT*1.0f );

    float xt = (xf*2.0 - 1.0f) * ( WIDTH*1.0f / HEIGHT);
    float yt = yf*2.0 - 1.0f;


    int samples = pImgDesc[0].sampleRate;

    for( int k = 0; k < samples; k++ ){

    float2 screen = findRandomPoint( &state, x*HEIGHT + y, (float2)( xt, yt ) );

    float4 ray = (shoot_ray( screen, pCamera ));

    FullIntersection patch;// = ray_intersect( pCamera[0].vPos, ray, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, -1, -1 );

    // Second intersection.


    FullIntersection eyePatch;
    eyePatch.vDir = ray;
    eyePatch.kInter.vPos = pCamera[0].vPos;
    eyePatch.uPrimType = -1;
    eyePatch.uIndex = -1;

    float4 col = (float4)(1.0f,1.0f,1.0f,1.0f);
    int i = 0;

    patch = eyePatch;

    float4 incomingDir = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

    for( ; i < 5; i++){


            float4 launchRef = patch.vDir;
            float4 launchPos = patch.kInter.vPos + patch.vDir * 0.00001f;// Push slightly to avoid self-intersection.
            //if( i == 0 )

            
            float4 launchDir;

    //         if( dot( launchRef, (launchPos - (float4)( 2.0f, 2.0f, 2.0f, 0.0f )) ) > 0  )
                launchDir = brdf_reflective( patch, (float4)(0,0,0,0), incomingDir, &state, x*HEIGHT + y );
    //         else
    //            launchDir = brdf_lambertian_imp( patch, (float4)(0,0,0,0), incomingDir, fast_normalize( -(launchPos - (float4)( 3.0f, 3.0f, 3.0f, 0.0f )) ), 0.0f, 0.2f * M_PI, &state, x*HEIGHT + y );

            float wt = launchDir.w;
            launchDir.w = 0.0f;

            // Directly launch the first ray( from the eye ).
            if( i == 0 ){
                wt = 1.0f;
                launchDir = launchRef;
            }

            if( i == 4 ){
                col *= (float4)(0.3f, 0.3f, 0.3f, 0.3f);// sky color.
                break;
            }

            //cmdlog("launchRef: %f, %f, %f, %f\n", launchRef.x, launchRef.y, launchRef.z, launchRef.w);
            //cmdlog("launchDir: %f, %f, %f, %f\n", launchDir.x, launchDir.y, launchDir.z, launchDir.w);
            
            FullIntersection patch2 = ray_intersect( launchPos, launchDir, pSpheres, pPlanes, pTriangles, pDesc, pSurfaces, patch.uPrimType, patch.uIndex );

            incomingDir = launchDir;

            if( patch2.kInter.bInter ){
                if( patch2.vEmissive.w != 0.0f ){
                    cmdlog( "weight: %f\n", wt );
                    col *= ((float4)patch2.vEmissive/wt) * ( (float)M_1_PI * 0.5f );// sky color.
                    cmdlog("%d HIT LIGHT %d,%d\n", i, patch2.uPrimType, patch2.uIndex);
                    break;
                }
                cmdlog( "weight: %f\n", wt );
                col *= ( (float4)patch2.vDiffuse/wt * (float)M_1_PI * 0.5f );
                cmdlog("%d HIT %d,%d\n", i, patch2.uPrimType, patch2.uIndex);
            }else{
                col *= (float4)(0.3f, 0.3f, 0.4f, 0.3f);// sky color.
                cmdlog("%d Miss\n", i);
                break;
            }
                //col = patch.vDiffuse * col;
        patch = patch2;
    }

    cmdlog("final col: %f, %f, %f, %f; %d bumps\n", col.x, col.y, col.z, col.w, i);

    float4 newSample = col;
    float4 currSample = prevValue * native_divide( ((float)numSamples), (numSamples + 1)) + newSample * native_divide( 1.0f, (float)(numSamples + 1) );
    if( numSamples == 0 )
      currSample = newSample;

    prevValue = currSample;
    numSamples ++;
    }

    randSeed[x*HEIGHT + y] = state;

    write_imagef( imgOut, (int2)(x,y), prevValue );
}*/
