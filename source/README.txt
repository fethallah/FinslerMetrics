
This directory contains the source files for the computing the multi-scale oriented flux matrix and some isotropic metrics associated to it. 

The folder data contains the sample images(tif, jpg or png).

The folder outputs is an empty folder which will hold the user's tests results.


In order to use this application you need to have the following software installed:
CMake 2.2.8 (or above)
Insight Toolkit 4.3, with FFTWD and FFTWF flags activated.

To run the provided synthetic example, you should build the project,
then, run the following line:
$> itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Synthetic-02.png 0 0 outputs/Synthetic-02_tubularity.nrrd 1 outputs/Synthetic-02_scale_space_tubularity.nrrd 0 0 0 0 0 0 4 8 9 0.5 1 1 100000
