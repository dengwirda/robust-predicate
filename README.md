##  `Robust geometric predicates (without the agonising pain)`

Evaluating "geometric predicates" (`orientation`, `incircle`, etc) using standard floating-point arithmetic is a well-known nightmare, with floating-point round-off errors leading to ambiguities, non-convergence, and program crashes. A standard remedy is to employ "exact" computation via multi-precision techniques; eliminating floating-point round-off errors through operations on arbitrary bit-length number types.

This package is a `C++` framework for the construction of such predicates; encapsulating <a href=https://www.cs.cmu.edu/~quake/robust.html>Jonathan Shewchuk's</a> seminal arbitrary precision library <a href=https://dl.acm.org/doi/book/10.5555/865018>*without the agonising pain*</a> of the original hand-rolled `C89` code.

This package aims to implement a "zero-overhead" abstraction; leveraging various `C++` template- and compile-time patterns to avoid run-time stack/heap manipulation or pointer indirection. Timing analysis suggests the new templated `C++` predicates perform favourably compared to the original hand-rolled `C89` implementation. 

The following predicates are currently available:
````
orient2d: orientation of 3 points in E^2, or a point wrt. a line.
bisect2d: orientation of point wrt. half-space in E^2.
bisect2w: orientation of point wrt. half-space in E^2 (weighted).
inball2d: point-in-circumball (Delaunay-Voronoi tessellations) in E^2. 
inball2w: point-in-ortho-ball (Regular-Laguerre tessellations) in E^2.

orient3d: orientation of 4 points in E^3, or a point wrt. a plane.
bisect3d: orientation of point wrt. half-space in E^3.
bisect3w: orientation of point wrt. half-space in E^3 (weighted).
inball3d: point-in-circumball (Delaunay-Voronoi tessellations) in E^3.
inball3w: point-in-ortho-ball (Regular-Laguerre tessellations) in E^3.
````
A simplified two-stage variation on <a href=https://doi.org/10.1007/PL00009321>Shewchuk's original arithmetic</a> is employed, adopting standard (fast!) floating-point approximations when results are unambiguous and falling back onto (slower) arbitrary precision evaluations as necessary to guarantee "sign-correctness". Semi-static filters are used to toggle between floating-point and arbitrary precision kernels.

### `License`

This program may be freely redistributed under the condition that the copyright notices (including this entire header) are not removed, and no compensation is received through use of the software.  Private, research, and institutional use is free.  You may distribute modified versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS`. Distribution of this code as part of a commercial system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR`. (If you are not directly supplying this code to a customer, and you are instead telling them how they can obtain it for free, then you are not required to make any arrangement with me.) 

`DISCLAIMER`:  Neither I nor: Columbia University, the Massachusetts Institute of Technology, the University of Sydney, nor the National Aeronautics and Space Administration warrant this code in any way whatsoever.  This code is provided "as-is" to be used at your own risk.

### `References`

`[1]` - J. R. Shewchuk (1997), Adaptive Precision Floating-Point Arithmetic & Fast Robust Geometric Predicates. Discrete & Computational Geometry, 18, pp. 305-363.

`[2]` - B. Lévy (2016), Robustness and efficiency of geometric programs: The Predicate Construction Kit (PCK). Computer-Aided Design, 72, pp. 03-12.

`[3]` - C. Burnikel, S. Funke, and M. Seel (2001), Exact geometric computation using cascading. IJCGA (Special issue) 11 (3), pp. 245–266.

