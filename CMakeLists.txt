cmake_minimum_required(VERSION 2.8)
project(bsplinefitting)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

add_definitions(-std=c++11)

# This looks for the libnabo-config.cmake file, which in turns
# includes libnabo-targets.cmake:
find_package(libnabo REQUIRED)

cs_add_library(${PROJECT_NAME}
  src/cubic_b_spline.cpp
  src/spline_curve_fitting.cpp
  src/open_cubic_b_spline.cpp
)
target_link_libraries(${PROJECT_NAME})

cs_install()
cs_export()

add_executable(bspline_fitting_example main.cpp )
target_link_libraries(bspline_fitting_example ${PROJECT_NAME} ${catkin_LIBRARIES} libnabo::nabo)