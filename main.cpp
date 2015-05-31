/*
  OpenCL entry point
*/

#include "CL/cl.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>


inline void checkErr(cl_int err, const char * name) {
    if (err != CL_SUCCESS) {
      std::cout << "ERROR: " << name  << " (" << err << ")" << std::endl;
      exit(EXIT_FAILURE);
   }
}

int main(void){
   cl_int err;

   std::vector<cl::Platform> platformList;
   cl::Platform::get(&platformList);

   checkErr(platformList.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");

   std::cerr << "Platform number is: " << platformList.size() << std::endl;
   std::string platformVendor;

   platformList[0].getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
   std::cerr << "Platform is by: " << platformVendor << "\n";

   cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platformList[0])(), 0};

   cl::Context context(
      CL_DEVICE_TYPE_CPU,
      cprops,
      NULL,
      NULL,
      &err);

   checkErr(err, "Context::Context()");

   std::cout<<"Created context."<< std::endl;

   const int inSize = 200000;

   float* outH = new float[inSize];
   cl::Buffer outCL( context, CL_MEM_WRITE_ONLY, inSize * sizeof( float ) );

   checkErr(err, "Buffer::Buffer()");

   cl::Buffer inCL( context, CL_MEM_READ_ONLY, inSize * sizeof( float ) );

   checkErr(err, "Buffer::Buffer()");


   std::cout<<"Created buffers."<< std::endl;

   float *f = new float[inSize];
   for( int i = 0; i < inSize; i++ ){
     f[i] = i;
   }


   std::vector<cl::Device> devices;
   devices = context.getInfo<CL_CONTEXT_DEVICES>();
   checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

   std::cout<<"Num available devices: "<< devices.size()<< std::endl;

   std::ifstream file("kernel0.cl");
   checkErr(file.is_open() ? CL_SUCCESS:-1, "kernel0.cl");
   std::string prog(
     std::istreambuf_iterator<char>(file),
     (std::istreambuf_iterator<char>()));

   cl::Program::Sources source(1,
     std::make_pair(prog.c_str(), prog.length()+1));

   cl::Program program(context, source);
   err = program.build(devices,"");

   checkErr(err, "Program::build()");
   std::cout<<"Built program"<< std::endl;


   cl::Kernel kernel(program, "square", &err);
   checkErr(err, "Kernel::Kernel()");
   err = kernel.setArg(0, outCL);
   checkErr(err, "Kernel::setArg()");
   err = kernel.setArg(1, inCL);
   checkErr(err, "Kernel::setArg()");


   std::cout<<"Built Kernel"<< std::endl;


   cl::CommandQueue queue(context, devices[0], 0, &err);
   checkErr(err, "CommandQueue::CommandQueue()");

   queue.enqueueWriteBuffer( inCL, CL_TRUE, 0, inSize * sizeof(float), f);

   cl::Event event;
   err = queue.enqueueNDRangeKernel(
     kernel,
     cl::NullRange,
     cl::NDRange(inSize),
     cl::NDRange(1, 1),
     NULL,
     &event);

   checkErr(err, "CommandQueue::enqueueNDRangeKernel()");
   std::cout<<"Kernel executing"<< std::endl ;

   event.wait();

   float *fout = new float[inSize];
   err = queue.enqueueReadBuffer(
      outCL,
      CL_TRUE,
      0,
      inSize * sizeof( float ),
      fout);

   checkErr(err, "ComamndQueue::enqueueReadBuffer()");
   std::cout << "Kernel output:" << std::endl;

   for( int i = 0; i < inSize; i++ ){
     std::cout<< fout[i] << std::endl;
   }

   std::cout<<"Kernel finished executing."<< std::endl;

   return EXIT_SUCCESS;
}
