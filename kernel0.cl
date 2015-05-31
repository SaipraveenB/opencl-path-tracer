#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable


__kernel void square(__global float * out, __global float * in) {
   size_t tid = get_global_id(0);

   out[tid] = in[tid] * in[tid];
}
