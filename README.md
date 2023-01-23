BSplineFitting
==============

Fitting cubic spline curve to 2d points

This is an implementation of the paper
["Fitting B-spline Curves to Point Clouds
by Curvature-Based Squared Distance Minimization"](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/12/Fitting-B-spline-Curves-to-Point-Clouds-by-Curvature-Based-Squared-Distance-Minimization.pdf) by Wang et al.

This fork puts the original repository into a catkin package and uses [libnabo](https://github.com/ethz-asl/libnabo) instead of [ANN](https://www.cs.umd.edu/~mount/ANN/).
```angular2html
# before catkin_make, you should git clone catkin_simple and eigen_catkin into your catkin workspace
# Please loop up the above two repositories on github
```
The package allows to fit open or closed B-spline curves through an unordered set of 2d points in contrast to the [curves package](https://github.com/ethz-asl/curves), which assumes points to be ordered in time.

```angular2html
#run
./bspline_fitting_example
```