/*
OpenCL entry point
*/

#define GLFW_EXPOSE_NATIVE_GLX
#define GLFW_EXPOSE_NATIVE_X11

#include <stdgl.h>

#ifdef __linux__
#include <GLFW/glfw3native.h>
#else
#include <OpenGL/OpenGL.h>
#endif

#include <assets/texture.h>
#include <renderer/render_target.h>
#include <renderer/camera.h>
#include <renderer/renderer.h>
#include <geometry/primitives.h>

#include <CL/cl.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <time.h>
#include <unistd.h>



#define WIDTH 1067
#define HEIGHT 600
#define SAMPLES 1

#define FPS_INTERVAL 0.5

inline void checkErr(cl_int err, const char * name) {
  if (err != CL_SUCCESS) {
    std::cout << "ERROR: " << name  << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}
inline void checkGLErr( const char* name ){
  GLint err = glGetError();
  if (err != GL_NO_ERROR) {
    std::cout << "ERROR: " << name  << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// GLFW window handle.
GLFWwindow* window;

// Main shared memory unit for rendering from OpenCL
RenderTarget* pCLTarget;

// Secondary shared memory to act as accumulator.
RenderTarget* pAccumulator;

// Vector of memory shared between OpenCL and OpenGL contexts.
std::vector<cl::Memory>* vSharedUnits;

// Camera class to control movement through the scene.
ModelCamera* pCamera;

ImageDescriptor imgDesc;

int iFrameCount = 0;
float fFrameAverage = 0.0f;
int iSamples = 0;

clock_t cLast;


void handleFrameCounter(){
  iFrameCount++;
  clock_t cNow = clock();

  //std::cout<<((float)(cNow - cLast)/(CLOCKS_PER_SEC))<<"\n";
  if( ((float)(cNow - cLast)/(CLOCKS_PER_SEC)) > FPS_INTERVAL ){
    cLast = cNow;
    fFrameAverage = (fFrameAverage * ( iSamples/(iSamples + 1)) + iFrameCount/( FPS_INTERVAL * (iSamples + 1) ) );
    iSamples++;

    std::cout<<"Differential: "<< ((float)iFrameCount)/FPS_INTERVAL<< std::endl;
    std::cout<<"Average: "<< fFrameAverage << " over " << iSamples << " samples." << std::endl;

    iFrameCount = 0;
  }

}

int sceneChanged(){
    return pCamera->isChanged();
}
void mainLoop( cl::CommandQueue& queue, cl::Context& context, cl::Kernel kernel, cl::Buffer clImgDesc, cl::Buffer clCamera ){
  cl::Event eAcquire, eRelease, eExecute;
  cl_int err;


  glFinish();
  checkGLErr( "glFinish()" );

  queue.enqueueWriteBuffer( clImgDesc, CL_TRUE, 0, 1 * sizeof(ImageDescriptor), (const void*)&imgDesc);

  err = queue.enqueueAcquireGLObjects( vSharedUnits, NULL, &eAcquire );
  checkErr(err, "CommandQueue::enqueueAcquireGLObjects()");

  eAcquire.wait();


  err = queue.enqueueNDRangeKernel( kernel, cl::NullRange, cl::NDRange(WIDTH, HEIGHT), cl::NullRange, NULL, &eExecute);

  checkErr(err, "CommandQueue::enqueueNDRangeKernel()");
  //std::cout<<"Kernel executing"<< std::endl ;
  clock_t ti = clock();
  eExecute.wait();
  clock_t tf = clock();

  queue.finish();
  err = queue.enqueueReleaseGLObjects( vSharedUnits, NULL, &eRelease );
  checkErr(err, "CommandQueue::enqueueReleaseGLObjects()");

  eRelease.wait();


  imgDesc.numSamples += SAMPLES;

  pAccumulator->glBind( GL_DRAW_FRAMEBUFFER );
  checkGLErr( "glBind GL_DRAW_FRAMEBUFFER, Accumulator " );
  pCLTarget->glBind( GL_READ_FRAMEBUFFER );
  checkGLErr( "glBind GL_READ_FRAMEBUFFER, Main Target " );
  glBlitFramebuffer( 0, 0, WIDTH, HEIGHT, 0, 0, WIDTH, HEIGHT, GL_COLOR_BUFFER_BIT, GL_NEAREST );
  checkGLErr( "glBlitFramebuffer" );

  glBindFramebuffer( GL_DRAW_FRAMEBUFFER, 0 );
  checkGLErr( "glBind GL_DRAW_FRAMEBUFFER, 0 " );
  pCLTarget->glBind( GL_READ_FRAMEBUFFER );
  checkGLErr( "glBind GL_READ_FRAMEBUFFER, something " );
  glBlitFramebuffer( 0, 0, WIDTH, HEIGHT, 0, 0, WIDTH, HEIGHT, GL_COLOR_BUFFER_BIT, GL_NEAREST );
  checkGLErr( "glBlitFramebuffer" );

  glfwPollEvents();

  pCamera->glfwHandleCursor( ((float)(tf - ti))/(CLOCKS_PER_SEC * 1.0f) );
  if( sceneChanged() ){
    //printf("scene changed..!");
    imgDesc.numSamples = 0;
    CLCamera* cam = pCamera->getCLCamera();
    queue.enqueueWriteBuffer( clCamera, CL_TRUE, 0, 1 * sizeof(CLCamera), (const void*)cam );
    delete cam;
  }

  glfwSwapBuffers( window );
  checkGLErr( "glSwapBuffers" );

  //Block for a while.
  //int i;
  //std::cin >> i;

  //float timeTaken = ( (float)(tf - ti) ) / (float)CLOCKS_PER_SEC;
  //std::cout<<"Time taken: "<< timeTaken * 1000 << "ms" << std::endl;
  //std::cout<<"Predicted FPS: "<< 1 / timeTaken << " FPS"<< std::endl;
  if( imgDesc.numSamples % 10 == 0 )
    std::cout<<"numSamples: "<<imgDesc.numSamples<<std::endl;
  //handleFrameCounter();

}



int main(int argc, char** argv){
  cl_int err;


  // START OPENGL INIT

  Magick::InitializeMagick(argv[0]);
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit ()) {
    fprintf (stderr, "ERROR: could not start GLFW3\n");
    return 1;
  }
  // Demand OpenGL 4.1
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

  #ifdef __APPLE__
  // Use Core profile to obtain a context for the latest OpenGL spec.
  // Otherwise we're stuck at 2.1
  std::cout<<"Apple FTW\n";
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  #endif




  window = glfwCreateWindow (WIDTH, HEIGHT, "Hello Triangle", NULL, NULL);
  if (!window) {
    fprintf (stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    return 1;
  }
  glfwSetWindowSize( window, WIDTH/2 , HEIGHT/2);

  //glfwWindowHint(GLFW_SAMPLES, 4);
  //glEnable( GL_MULTISAMPLE );

  glfwMakeContextCurrent (window);
  checkGLErr( "glfwMakeContextCurrent" );
  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  glewInit();
  glGetError();
  //checkGLErr( "GLEW init" );



  // END OPENGL INIT..

  // START OPENCL..
  std::vector<cl::Platform> platformList;
  cl::Platform::get(&platformList);

  checkErr(platformList.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");

  std::cerr << "Platform number is: " << platformList.size() << std::endl;
  std::string platformVendor;

  platformList[0].getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
  std::cerr << "Platform is by: " << platformVendor << "\n";

  #ifdef __APPLE__
  CGLContextObj kCGLContext = CGLGetCurrentContext();
  CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);
  cl_context_properties cprops[6] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platformList[0])(),CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE , (cl_context_properties) kCGLShareGroup, 0, 0};
  cl::Context context( CL_DEVICE_TYPE_CPU, cprops, NULL, NULL, &err);
  #endif

  #ifdef __linux__
    cl_platform_id platform;
    err = clGetPlatformIDs(1, &platform, NULL);
    cl_context_properties props[] =
    {
    	CL_GL_CONTEXT_KHR, (cl_context_properties)glfwGetGLXContext( window ),
    	CL_GLX_DISPLAY_KHR, (cl_context_properties)glfwGetX11Display(),
    	CL_CONTEXT_PLATFORM, (cl_context_properties)platform,
    	0
    };
    cl::Context context( CL_DEVICE_TYPE_CPU, props, NULL, NULL, &err);
  //cl::Context context = clCreateContextFromType(props, CL_DEVICE_TYPE_CPU, NULL, NULL, &err);
  #endif


  checkErr(err, "Context::Context()");

  std::cout<<"Created context."<< std::endl;

  // Create shared texture.
  pCLTarget = new RenderTarget( WIDTH, HEIGHT, GL_RGBA , GL_RGBA, GL_FLOAT, 0, false );
  checkGLErr( "RenderTarget::RenderTarget" );
  pAccumulator = new RenderTarget( WIDTH, HEIGHT, GL_RGBA, GL_RGBA, GL_FLOAT, 0, false );
  checkGLErr( "RenderTarget::RenderTarget" );

  const int inSizeS = 3;
  const int inSizeP = 0;
  const int inSizeT = 12;
  const int inSurf = 5;

  /*
  float* outH = new float[inSize];
  cl::Buffer outCL( context, CL_MEM_WRITE_ONLY, inSize * sizeof( float ) );
  */
  Sphere* spheres = new Sphere[inSizeS];
  std::cout<<"Sphere: "<< spheres[0].radius << "\n";
  spheres[1].uSurf = 0;
  spheres[1].center = glm::vec4( 0.0f, -2.0f, 0.0f, 0.0f );
  spheres[1].radius = 1.0f;
  /*spheres[1].uSurf = 0;
  spheres[1].center = glm::vec4( +1.0f, 0.0f, +1.0f, 0.0f);
  spheres[1].radius = 1.0f;*/

  spheres[0].uSurf = 2;
  spheres[0].center = glm::vec4( 0.0f, +1.50f, 0.0f, 0.0f);
  spheres[0].radius = 0.2f;

  /*spheres[2].uSurf = 0;
  spheres[2].center = glm::vec4( +0.0f, 0.0f, -0.0f, 0.0f);
  spheres[2].radius = 1.0f;*/

  spheres[2].uSurf = 0;
  spheres[2].center = glm::vec4( -2.0f, -2.0f, -2.0f, 0.0f);
  spheres[2].radius = 1.0f;

  /*spheres[4].uSurf = 0;
  spheres[4].center = glm::vec4( -1.0f, -0.0f, +1.0f, 0.0f);
  spheres[4].radius = 1.0f;*/
  Plane* planes = new Plane[inSizeP];
  int planeSize = 3.0f;
  //std::cout<<"Sphere: "<< planes[0].radius << "\n";
  /*planes[0].normal = glm::vec4( 0.0f, 1.0f, 0.0f, 0.0f );
  planes[0].point = glm::vec4( 0.0f, -planeSize, 0.0f, 0.0f );
  planes[0].uSurf = 1;
  
  planes[1].normal = glm::vec4( 0.0f, -1.0f, 0.0f, 0.0f );
  planes[1].point = glm::vec4( 0.0f, planeSize, 0.0f, 0.0f );
  planes[1].uSurf = 1;
  
  planes[2].normal = glm::vec4( 1.0f, 0.0f, 0.0f, 0.0f );
  planes[2].point = glm::vec4( -planeSize, 0.0f, 0.0f, 0.0f );
  planes[2].uSurf = 1;
  
  planes[3].normal = glm::vec4( -1.0f, 0.0f, 0.0f, 0.0f );
  planes[3].point = glm::vec4( planeSize, 0.0f, 0.0f, 0.0f );
  planes[3].uSurf = 1;
  
  planes[4].normal = glm::vec4( 0.0f, 0.0f, +1.0f, +0.0f );
  planes[4].point = glm::vec4( 0.0f, 0.0f, -planeSize, 0.0f );
  planes[4].uSurf = 1;
  
  planes[5].normal = glm::vec4( 0.0f, 0.0f, -1.0f, 0.0f );
  planes[5].point = glm::vec4( 0.0f, 0.0f, +planeSize, 0.0f );
  planes[5].uSurf = 1;*/
  
  

  Triangle* triangles = new Triangle[inSizeT];
  //std::cout<<"Sphere: "<< spheres[0].radius << "\n";
  float boxSize = 3.0f;
  triangles[0].uSurf = 1;
  triangles[0].p0 = glm::vec4( -boxSize, -boxSize, -boxSize, 0.0f);
  triangles[0].p1 = glm::vec4( -boxSize, boxSize, -boxSize, 0.0f);
  triangles[0].p2 = glm::vec4( -boxSize, -boxSize, boxSize, 0.0f);
  
  triangles[1].uSurf = 1;
  triangles[1].p0 = glm::vec4( -boxSize, boxSize, boxSize, 0.0f);
  triangles[1].p2 = glm::vec4( -boxSize, boxSize, -boxSize, 0.0f);
  triangles[1].p1 = glm::vec4( -boxSize, -boxSize, boxSize, 0.0f);
  
  triangles[2].uSurf = 1;
  triangles[2].p0 = glm::vec4( boxSize, -boxSize, -boxSize, 0.0f);
  triangles[2].p2 = glm::vec4( boxSize, boxSize, -boxSize, 0.0f);
  triangles[2].p1 = glm::vec4( boxSize, -boxSize, boxSize, 0.0f);
  
  triangles[3].uSurf = 1;
  triangles[3].p0 = glm::vec4( boxSize, boxSize, boxSize, 0.0f);
  triangles[3].p1 = glm::vec4( boxSize, boxSize, -boxSize, 0.0f);
  triangles[3].p2 = glm::vec4( boxSize, -boxSize, boxSize, 0.0f);
  

  triangles[4].uSurf = 1;
  triangles[4].p0 = glm::vec4( -boxSize, -boxSize, -boxSize, 0.0f);
  triangles[4].p2 = glm::vec4( boxSize, -boxSize, -boxSize, 0.0f);
  triangles[4].p1 = glm::vec4( -boxSize, -boxSize, boxSize, 0.0f);
  
  triangles[5].uSurf = 1;
  triangles[5].p0 = glm::vec4( boxSize, -boxSize, boxSize, 0.0f);
  triangles[5].p1 = glm::vec4( boxSize, -boxSize, -boxSize, 0.0f);
  triangles[5].p2 = glm::vec4( -boxSize, -boxSize, boxSize, 0.0f);
  
  triangles[6].uSurf = 1;
  triangles[6].p0 = glm::vec4( -boxSize, boxSize, -boxSize, 0.0f);
  triangles[6].p1 = glm::vec4( boxSize, boxSize, -boxSize, 0.0f);
  triangles[6].p2 = glm::vec4( -boxSize, boxSize, boxSize, 0.0f);
  
  triangles[7].uSurf = 1;
  triangles[7].p0 = glm::vec4( boxSize, boxSize, boxSize, 0.0f);
  triangles[7].p2 = glm::vec4( boxSize, boxSize, -boxSize, 0.0f);
  triangles[7].p1 = glm::vec4( -boxSize, boxSize, boxSize, 0.0f);
  
  triangles[8].uSurf = 3;
  triangles[8].p0 = glm::vec4( -boxSize, -boxSize, -boxSize, 0.0f);
  triangles[8].p1 = glm::vec4( boxSize, -boxSize, -boxSize, 0.0f);
  triangles[8].p2 = glm::vec4( -boxSize, boxSize, -boxSize, 0.0f);
  
  triangles[9].uSurf = 3;
  triangles[9].p0 = glm::vec4( boxSize, boxSize, -boxSize, 0.0f);
  triangles[9].p2 = glm::vec4( boxSize, -boxSize, -boxSize, 0.0f);
  triangles[9].p1 = glm::vec4( -boxSize, boxSize, -boxSize, 0.0f);
  
  triangles[10].uSurf = 4;
  triangles[10].p0 = glm::vec4( -boxSize, -boxSize, boxSize, 0.0f);
  triangles[10].p2 = glm::vec4( boxSize, -boxSize, boxSize, 0.0f);
  triangles[10].p1 = glm::vec4( -boxSize, boxSize, boxSize, 0.0f);
  
  triangles[11].uSurf = 4;
  triangles[11].p0 = glm::vec4( boxSize, boxSize, boxSize, 0.0f);
  triangles[11].p1 = glm::vec4( boxSize, -boxSize, boxSize, 0.0f);
  triangles[11].p2 = glm::vec4( -boxSize, boxSize, boxSize, 0.0f);
  
  /*triangles[12].uSurf = 1;
  triangles[12].p0 = glm::vec4( -boxSize/4.0f, boxSize, -boxSize/4.0f, 0.0f);
  triangles[12].p1 = glm::vec4( boxSize/4.0f, boxSize, -boxSize/4.0f, 0.0f);
  triangles[12].p2 = glm::vec4( -boxSize/4.0f, boxSize, boxSize/4.0f, 0.0f);
  
  triangles[13].uSurf = 1;
  triangles[13].p0 = glm::vec4( boxSize/4.0f, boxSize, boxSize/4.0f, 0.0f);
  triangles[13].p2 = glm::vec4( boxSize/4.0f, boxSize, -boxSize/4.0f, 0.0f);
  triangles[13].p1 = glm::vec4( -boxSize/4.0f, boxSize, boxSize/4.0f, 0.0f);*/
  
  
  

  GeometryDescriptor* geometry = new GeometryDescriptor( inSizeS, inSizeP, inSizeT );

  Surface* pSurf = new Surface[inSurf];
  pSurf[0].vColor = glm::vec4( 1.0f, 1.0f, 0.0f, 1.0f );
  
  pSurf[1].vColor = glm::vec4( 1.0f, 1.0f, 1.0f, 1.0f );
  
  pSurf[2].vColor = glm::vec4( 1.0f, 1.0f, 1.0f, 1.0f );
  float lPower = 5.0f;
  pSurf[2].vEmissive = glm::vec4( lPower, lPower, lPower, lPower );
  
  pSurf[3].vColor = glm::vec4( 1.0f, 0.0f, 0.0f, 1.0f );
  
  pSurf[4].vColor = glm::vec4( 0.0f, 1.0f, 0.0f, 1.0f );

  cl::Buffer clSpheres( context, CL_MEM_READ_ONLY, inSizeS * sizeof( Sphere ));
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clPlanes( context, CL_MEM_READ_ONLY, inSizeP * sizeof( Plane ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clTriangles( context, CL_MEM_READ_ONLY, inSizeT * sizeof( Triangle ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clCamera( context, CL_MEM_READ_ONLY, 1 * sizeof( CLCamera ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clGeom( context, CL_MEM_READ_ONLY, 1 * sizeof( GeometryDescriptor ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clImgDesc( context, CL_MEM_READ_ONLY, 1 * sizeof( ImageDescriptor ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clSeed( context, CL_MEM_READ_WRITE, WIDTH * HEIGHT * 4 * sizeof( uint ) );
  checkErr(err, "Buffer::Buffer()");

  cl::Buffer clSurf( context, CL_MEM_READ_WRITE, inSurf * sizeof( Surface ) );
  checkErr(err, "Buffer::Buffer()");

  cl::ImageGL imgGL( context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, pCLTarget->getColorTexture()->glGetInternalTexture(), &err );
  checkErr(err, "ImageGL::ImageGL()");

  cl::ImageGL accGL( context, CL_MEM_READ_ONLY, GL_TEXTURE_2D, 0, pAccumulator->getColorTexture()->glGetInternalTexture(), &err );
  checkErr(err, "ImageGL::ImageGL()");

  std::cout<<"Created buffers."<< std::endl;

  srand( time( NULL ) );
  uint *pSeeds = new uint[WIDTH * HEIGHT * 4];
  for( int i = 0; i < WIDTH * HEIGHT * 4; i++ ){
    pSeeds[i] = rand();
  }

  std::vector<cl::Device> devices;
  devices = context.getInfo<CL_CONTEXT_DEVICES>();
  checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

  std::cout<<"Num available devices: "<< devices.size()<< std::endl;

  std::ifstream file("src/kernel/kernel0.cl");
  checkErr(file.is_open() ? CL_SUCCESS:-1, "src/kernel/kernel0.cl");
  std::string prog( std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));

  cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length()+1));

  std::cout<<"Source obtained."<< std::endl;
  
  cl::Program program(context, source);
  err = program.build(devices,"-cl-opt-disable");
  std::cout<<"Source obtained."<< std::endl;

  std::string buildLog;
  program.getBuildInfo( devices[0], CL_PROGRAM_BUILD_LOG, &buildLog );
  std::cout<<"Build log:" << buildLog<< std::endl;
  checkErr(err, "Program::build()");
  std::cout<<"Built program"<< std::endl;

  cl::Kernel kernel(program, "bi_directional_path_trace", &err);
  checkErr(err, "Kernel::Kernel()");

  err = kernel.setArg(0, clCamera);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(1, imgGL);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(2, accGL);
  checkErr(err, "Kernel::setArg()");

  err = kernel.setArg(3, clSpheres);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(4, clPlanes);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(5, clTriangles);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(6, clGeom);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(7, clImgDesc);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(8, clSeed);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(9, clSurf);
  checkErr(err, "Kernel::setArg()");

  std::cout<<"Built Kernel"<< std::endl;

  pCamera = new ModelCamera( window );
  pCamera->setSpeedX( 0.03f );
  pCamera->setSpeedY( 0.03f );

  pCamera->setRadius( 8.0f );
  pCamera->setOrientation( glm::vec3( 0.0f, -1.0f, 0.0f ) );
  pCamera->reset( glm::vec3( 1.0f, 0.1f, -0.1f ) );


  cl::CommandQueue queue(context, devices[0], 0, &err);
  checkErr(err, "CommandQueue::CommandQueue()");
  CLCamera* cam = pCamera->getCLCamera();

  std::cout<<cam->vPos.x<<","<<cam->vPos.y<<","<<cam->vPos.z<<std::endl;

  std::cout<<cam->vLookAt.x<<","<<cam->vLookAt.y<<","<<cam->vLookAt.z<<std::endl;

  std::cout<<cam->vUp.x<<","<<cam->vUp.y<<","<<cam->vUp.z<<std::endl;

  std::cout<< sizeof( Plane )<< std::endl;
  queue.enqueueWriteBuffer( clCamera, CL_TRUE, 0, 1 * sizeof(CLCamera), (const void*)cam );
  queue.enqueueWriteBuffer( clSpheres, CL_TRUE, 0, inSizeS * sizeof(Sphere), (const void*)spheres);
  queue.enqueueWriteBuffer( clPlanes, CL_TRUE, 0, inSizeP * sizeof(Plane), (const void*)planes);
  queue.enqueueWriteBuffer( clTriangles, CL_TRUE, 0, inSizeT * sizeof(Triangle), (const void*)triangles);
  queue.enqueueWriteBuffer( clGeom, CL_TRUE, 0, 1 * sizeof(GeometryDescriptor), (const void*)geometry);
  queue.enqueueWriteBuffer( clSeed, CL_TRUE, 0, WIDTH * HEIGHT * 4 * sizeof(uint), (const void*)pSeeds);
  queue.enqueueWriteBuffer( clSurf, CL_TRUE, 0, inSurf * sizeof(Surface), (const void*)pSurf);

  vSharedUnits = new std::vector<cl::Memory>();
  vSharedUnits->push_back( imgGL );
  vSharedUnits->push_back( accGL );

  //Initialise counter.
  imgDesc.numSamples = 0;
  imgDesc.sampleRate = SAMPLES;
  cLast = clock();
  while( !glfwWindowShouldClose( window ) ){
    //usleep( 1000000 );
    mainLoop( queue, context, kernel, clImgDesc, clCamera );
  }

  /* Previous Program. Remove these if you think they are not required.
  float *fout = new float[inSize];
  err = queue.enqueueReadBuffer( clSpheres, CL_TRUE, 0, inSize * sizeof(Sphere), fout);
  */
  checkErr(err, "ComamndQueue::enqueueReadBuffer()");

  std::cout<<"Kernel finished executing."<< std::endl;

  delete vSharedUnits;
  return EXIT_SUCCESS;
}
