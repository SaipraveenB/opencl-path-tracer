# opencl-path-tracer

This GPU path tracer implements the Bi Directional Path Tracing algorithm.

It also lays out the basic framework for all tracer type algorithms( convenience functions for shifting resources between the CPU & GPU. )

The algorithm works with simple primitives. Compared to traditional rasterizer graphics systems, tracer graphics support a much wider range of primitives without much computational overhead( anything with a closed form equation ).

## Primitives:

Current Primitives:
- Sphere
- Plane( Infinite )
- Triangle

Primitives that can be implemented within the current framework:
( All shapes with a fast ray-intersection test )
- Ellipsoid( Consequently, Spheres too )
- Elliptic cylinder( Consequently, Cylinders too ) 
- Plane( Infinite & Axis-aligned truncations )
- Cuboids
- Triangles

The implementation is interactive( not usable in real-time though as it still takes too long to converge ). It incrementally improves the image quality when there is no camera movement; 

Refreshes the entire buffer when the camera is moved.

Rendering of the cornell box with two yellow spheres(after about a 2 seconds )
![I_1 NEAR PATH TRACING_2.png](/I_1 NEAR PATH TRACING_2.png)

( The black spot is the light source ).

Rendering of the cornell box with two yellow spheres(after about 8 seconds )
![I_1 NEAR PATH TRACING.png](/I_1 NEAR PATH TRACING.png)

## Numerical issues.
Implementing the original path tracer algorithm with inverse square illumination while computing the path weights in BDPT leads to bright spots called 'fireflies' which occur when a contiguous region has surface patches with all the illumination coming from a disproportionately small part of the domain.

This makes convergence a lot slower since a small number of points in the region hit the jackpot and get lit up whereas the majority remains dark giving the false impression of randomly lit spots.

An example of an extreme case:
![Screen Shot 2015-06-16 at 10.30.22 pm.png](/Screen Shot 2015-06-16 at 10.30.22 pm.png)

A slightly more modest case( it's still irritating and not an acceptable rendering ).

The fireflies in the image were inflenced by a collection of numeric problems of which the inverse square was the biggest problem.
![Screen Shot 2015-06-16 at 3.50.51 pm.png](/Screen Shot 2015-06-16 at 3.50.51 pm.png)


A laughably simple fix is to simply bound the maximum lighting from any single connection. This introduces bias into the lighting estimator but still provides  acceptable global illumination.

An image of indirect illumination only:
![Screen Shot 2015-06-16 at 10.40.57 pm.png](/Screen Shot 2015-06-16 at 10.40.57 pm.png)

After some other simple numerical fixes, the image turns out like this:
![Screen Shot 2015-06-16 at 3.45.25 pm.png](/Screen Shot 2015-06-16 at 3.45.25 pm.png)

Several other numeric issues arise that have been tackled to give a speedy convergence to a beleivable image( retaining the global illumination property )

The accepted method to reduce 'fireflies' and increase convergence speed is Importance Sampling as it acheives the same effect without introducing bias.

## Acceleration Structures
Currently, the implementation for speeding up geometry intersections through space partitioning is somewhat naive.
It uses some simple math for culling objects behind the ray origin, but otherwise an exhaustive search of the object space is performed.

An object-space Octree plus bounding boxes provide considerable acceleration for incredibly complex scenes. The implementation is underway for both building and using an Octree using OpenCL to take advantage of any parallelizable points.

## Rendering
Outdoor rendering to observe light bleeding out of the holes in spherical arrangement.
![Screen Shot 2015-06-11 at 6.20.36 pm.png](/Screen Shot 2015-06-11 at 6.20.36 pm.png)
![Screen Shot 2015-06-11 at 6.21.07 pm copy.png](/Screen Shot 2015-06-11 at 6.21.07 pm copy.png)

Indoor illumination: Almost completely converged image.( The red hue is due to reflections off the adjacent red wall. ). The image was rendered with exaggerated anti-alisaing which is why edges seem almost non-existent but the image is slightly blurry.
![Screen Shot 2015-06-14 at 12.18.15 am copy.png](/Screen Shot 2015-06-14 at 12.18.15 am copy.png)
