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

![]({{site.baseurl}}//I_1%20NEAR%20PATH%20TRACING_2.png)
( The black spot is the light source ).

Rendering of the cornell box with two yellow spheres(after about 8 seconds )

![]({{site.baseurl}}//I_1%20NEAR%20PATH%20TRACING.png)

## Numerical issues.
Implementing the original path tracer algorithm with inverse square illumination while computing the path weights in BDPT leads to bright spots called 'fireflies' which occur when a contiguous region has surface patches with all the illumination concentrated in a narrow solid angle.

This makes convergence a lot slower since a small number of points in the region hit the jackpot and get lit up whereas the majority remains dark giving the false impression of randomly lit spots.

A laughably simple fix is to simply bound the maximum lighting from any connection. This affects the unbiased nature of lighting but still provides  global illumination.

Unfortunately, this form of a truncated inverse square function gives visible bands because while human perception can't tell apart similar-looking colours when viewed separately, it can perceive a sharp change in gradient.

Another simple fix is to smooth the function with something like this:
![]({{site.baseurl}}/INV_SQ%20VS%20MODIFIED_INV_SQ_2.png)

