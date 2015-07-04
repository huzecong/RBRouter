cmake_minimum_required(VERSION 3.0.2)
project(RBRouter)

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
	if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
		cmake_policy(VERSION 2.8.4)
	else()
		cmake_policy(VERSION 2.6)
	endif()
endif()

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif()

# CGAL and its components
find_package(CGAL QUIET COMPONENTS)

if(NOT CGAL_FOUND)
	message(STATUS "This project requires the CGAL library, and will not be compiled.")
	return()
endif()

# include helper file
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)
	message(STATUS "This project requires the Boost library, and will not be compiled.")
	return()
endif()

# include for local directory

# include for local package

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")

set(SOURCE_FILES src/main.cpp src/RBRouter.cpp src/RBSequentialEmbedding.cpp
src/RBBase.cpp src/RBShortestPath.cpp src/Box2D/b2BlockAllocator.cpp
src/Box2D/b2TrackedBlock.cpp src/Box2D/b2Settings.cpp src/HYPlotPS/HYPlotPS.cpp)
add_executable(RBRouter ${SOURCE_FILES})