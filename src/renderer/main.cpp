/*
OpenCL entry point
*/
#include <stdgl.h>
#include <OpenGL/OpenGL.h>
#include <assets/texture.h>
#include <renderer/render_target.h>

#include <CL/cl.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <time.h>

#include "sphere.h"

#define WIDTH 1920
#define HEIGHT 1080
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

// Vector of memory shared between OpenCL and OpenGL contexts.
std::vector<cl::Memory>* vSharedUnits;


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
void mainLoop( cl::CommandQueue& queue, cl::Context& context, cl::Kernel kernel ){
  cl::Event eAcquire, eRelease, eExecute;
  cl_int err;


  glFinish();
  checkGLErr( "glFinish()" );

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

  glBindFramebuffer( GL_DRAW_FRAMEBUFFER, 0 );
  checkGLErr( "glBind GL_DRAW_FRAMEBUFFER, 0 " );
  pCLTarget->glBind( GL_READ_FRAMEBUFFER );
  checkGLErr( "glBind GL_READ_FRAMEBUFFER, something " );

  glBlitFramebuffer( 0, 0, WIDTH, HEIGHT, 0, 0, WIDTH, HEIGHT, GL_COLOR_BUFFER_BIT, GL_NEAREST );
  checkGLErr( "glBlitFramebuffer" );

  glfwPollEvents();

  glfwSwapBuffers( window );
  checkGLErr( "glSwapBuffers" );

  //Block for a while.
  //int i;
  //std::cin >> i;

  //float timeTaken = ( (float)(tf - ti) ) / (float)CLOCKS_PER_SEC;
  //std::cout<<"Time taken: "<< timeTaken * 1000 << "ms" << std::endl;
  //std::cout<<"Predicted FPS: "<< 1 / timeTaken << " FPS"<< std::endl;

  handleFrameCounter();

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
  glewExperimental = TRUE;
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
  #endif

  cl_context_properties cprops[6] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platformList[0])(),CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE , (cl_context_properties) kCGLShareGroup, 0, 0};

  cl::Context context( CL_DEVICE_TYPE_CPU, cprops, NULL, NULL, &err);

  checkErr(err, "Context::Context()");

  std::cout<<"Created context."<< std::endl;

  // Create shared texture.
  pCLTarget = new RenderTarget( WIDTH, HEIGHT, GL_RGBA, GL_RGBA, GL_FLOAT, 0, false );
  checkGLErr( "RenderTarget::RenderTarget" );

  const int inSize = 2000000;
  /*
  float* outH = new float[inSize];
  cl::Buffer outCL( context, CL_MEM_WRITE_ONLY, inSize * sizeof( float ) );
  */
  Sphere* spheres = new Sphere[inSize];
  cl::Buffer clSpheres( context, CL_MEM_WRITE_ONLY, inSize * sizeof(Sphere));

  checkErr(err, "Buffer::Buffer()");
  cl::Buffer clCamera( context, CL_MEM_READ_ONLY, inSize * sizeof( float ) );

  checkErr(err, "Buffer::Buffer()");

  cl::ImageGL imgGL( context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, pCLTarget->getColorTexture()->glGetInternalTexture(), &err );

  checkErr(err, "ImageGL::ImageGL()");

  std::cout<<"Created buffers."<< std::endl;

  float *f = new float[inSize];
  for( int i = 0; i < inSize; i++ ){
    f[i] = i;
  }

  std::vector<cl::Device> devices;
  devices = context.getInfo<CL_CONTEXT_DEVICES>();
  checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

  std::cout<<"Num available devices: "<< devices.size()<< std::endl;

  std::ifstream file("src/kernel/kernel0.cl");
  checkErr(file.is_open() ? CL_SUCCESS:-1, "src/kernel/kernel0.cl");
  std::string prog( std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));

  cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length()+1));

  cl::Program program(context, source);
  err = program.build(devices,"");

  std::string buildLog;
  program.getBuildInfo( devices[0], CL_PROGRAM_BUILD_LOG, &buildLog );
  std::cout<<"Build log:" << buildLog<< std::endl;
  checkErr(err, "Program::build()");
  std::cout<<"Built program"<< std::endl;

  cl::Kernel kernel(program, "square", &err);
  checkErr(err, "Kernel::Kernel()");

  err = kernel.setArg(0, clSpheres);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(1, clCamera);
  checkErr(err, "Kernel::setArg()");
  err = kernel.setArg(2, imgGL);
  checkErr(err, "Kernel::setArg()");


  std::cout<<"Built Kernel"<< std::endl;


  cl::CommandQueue queue(context, devices[0], 0, &err);
  checkErr(err, "CommandQueue::CommandQueue()");

  queue.enqueueWriteBuffer( clCamera, CL_TRUE, 0, inSize * sizeof(float), f);
  queue.enqueueWriteBuffer( clSpheres, CL_TRUE, 0, inSize * sizeof(Sphere), f);

  vSharedUnits = new std::vector<cl::Memory>();
  vSharedUnits->push_back( imgGL );

  //Initialise counter.

  cLast = clock();
  while( !glfwWindowShouldClose( window ) ){
    mainLoop( queue, context, kernel );
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
