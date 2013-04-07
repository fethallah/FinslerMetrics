cmake_minimum_required(VERSION 2.8)
#Change PROJECT_NAME to the name of your project
PROJECT(FinslerianMetrics)

#include some macros from another file...
INCLUDE(IJMacros.txt)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

find_package(ITK 4.1 REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
endif()

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
endif()


set(CMAKE_BUILD_TYPE release) #debug none release distribution

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)


#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
${ITK_LIBRARIES} 
)
INCLUDE(${ITK_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR})

SET(CurrentExe "itkMultiScaleOrientedFluxBasedMeasureImageFilter")
ADD_EXECUTABLE(itkMultiScaleOrientedFluxBasedMeasureImageFilter itkMultiScaleOrientedFluxBasedMeasureImageFilter.cxx)
TARGET_LINK_LIBRARIES(itkMultiScaleOrientedFluxBasedMeasureImageFilter ${Libraries} boost_math_tr1)

########################################## testing ##########################################
####### Synthetic-02 #######
ADD_TEST(tubularity2D-Synthetic-02 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Synthetic-02.png 0 0 outputs/Synthetic-02_tubularity.nrrd 1 outputs/Synthetic-02_scale_space_tubularity.nrrd 0 0 0 0 0 0 4 8 9 0.5 1 1 100000
)

####### Real-DSA-01 #######

ADD_TEST(tubularity2D-Real-DSA-01 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Real-DSA-01.jpg 0 0 outputs/Real-DSA-01_tubularity.nrrd 1 outputs/Real-DSA-01_scale_space_tubularity.nrrd 0 0 0 0 0 0 1 9 17 0.5 0 1 100000)

####### OP_1 #######

ADD_TEST(tubularity3D-OP_1 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/OP_1.tif 0 0 outputs/OP_1_tubularity.nrrd 1 outputs/OP_1_scale_space_tubularity.nrrd 0 0 0 0 0 0 1 6 11 0.5 1 1 100000)

########################################## testing done ##########################################












