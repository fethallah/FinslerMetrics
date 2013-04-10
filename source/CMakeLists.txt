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

SET(CurrentExe "GenerateRiemannianMetrics")
ADD_EXECUTABLE(GenerateRiemannianMetrics GenerateRiemannianMetrics.cxx)
TARGET_LINK_LIBRARIES(GenerateRiemannianMetrics ${Libraries})

########################################## testing ##########################################
####### Synthetic-02 #######
ADD_TEST(tubularity2D-Synthetic-02 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Synthetic-02.png outputs/Synthetic-02_tubularity.nii 1 outputs/Synthetic-02_scale_space_tubularity.nii 1 outputs/Synthetic-02_OF_Matrix.nii 1 outputs/Synthetic-02_scale_space_OF_Matrix.nii 0 0 4 8 9 0.5 1 1 100000
)

ADD_TEST(Riemannian-Metric-Space-Synthetic-02 GenerateRiemannianMetrics outputs/Synthetic-02_OF_Matrix.nii outputs/Synthetic-02_Riemannian_Metric.nii 1 10 1 1 outputs/Synthetic-02_Tensor_Glyph.vtk 0.005)


####### Real-DSA-01 #######

ADD_TEST(tubularity2D-Real-DSA-01 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Real-DSA-01.jpg outputs/Real-DSA-01_tubularity.nii 1 outputs/Real-DSA-01_scale_space_tubularity.nii 1 outputs/Real-DSA-01_OF_Matrix.nii 1 outputs/Real-DSA-01_scale_space_OF_Matrix.nii 0 0 1 9 17 0.5 0 1 100000)

ADD_TEST(Riemannian2D-Real-DSA-01 GenerateRiemannianMetrics outputs/Real-DSA-01_OF_Matrix.nii outputs/Real-DSA-01_Riemannian_Metric.nii 0 10 1 1 outputs/Real-DSA-01__Tensor_Glyph.vtk 0.005)

####### Real-DSA-02 #######

ADD_TEST(tubularity2D-Real-DSA-02 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Real-DSA-02.tif outputs/Real-DSA-02_tubularity.nii 1 outputs/Real-DSA-02_scale_space_tubularity.nii 1 outputs/Real-DSA-02_OF_Matrix.nii 1 outputs/Real-DSA-02_scale_space_OF_Matrix.nii 0 0 1 9 17 0.5 0 1 100000)

ADD_TEST(Riemannian2D-Real-DSA-02 GenerateRiemannianMetrics outputs/Real-DSA-02_OF_Matrix.nii outputs/Real-DSA-02_Riemannian_Metric.nii 0 10 1 1 outputs/Real-DSA-02__Tensor_Glyph.vtk 0.05)

####### OP_1 #######

ADD_TEST(tubularity3D-OP_1 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/OP_1.tif outputs/OP_1_tubularity.nii 1 outputs/OP_1_scale_space_tubularity.nii 1 outputs/OP_1_OF_Matrix.nii 1 outputs/OP_1_scale_space_OF_Matrix.nii 0 0 1 6 11 0.5 1 1 100000)

ADD_TEST(Riemannian3D-OP_1 GenerateRiemannianMetrics outputs/OP_1_OF_Matrix.nii outputs/OP_1_OF_Riemannian_Metric.nii 1 2 1 1 outputs/OP_1_TensorGlyphs.vtk 0.001)

########################################## testing done ##########################################












