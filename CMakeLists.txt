cmake_minimum_required(VERSION 3.1.2)
project(mas_factor)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  geometry_msgs
  )

find_package(GTSAM REQUIRED)

set(LIBRARIES
  MasFactor
  )

catkin_package(
  INCLUDE_DIRS include
  DEPENDS GTSAM

  CATKIN_DEPENDS
  geometry_msgs

  LIBRARIES ${LIBRARIES}
  )

# add definitions
add_definitions(
  ${GTSAM_DEFINITIONS}
  )

# include directories
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${GTSAM_INCLUDE_DIRS}
  )

# link directories
link_directories(
  include
  ${GTSAM_LIBRARY_DIRS}
  )

###########
## Build ##
###########

# Mas Factor
add_library(MasFactor src/mas_factor.cpp src/mas_preintegration.cpp src/mas_preintegration_base.cpp src/full_state.cpp src/mas_bias.cpp)
add_dependencies(MasFactor
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
  )

target_link_libraries(MasFactor
  ${catkin_LIBRARIES}
  ${GTSAM_LIBRARIES}
  gtsam
  )

