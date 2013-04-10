//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// Copyright (C) 2013  Fethallah Benmansour                                     //
//																				//
// This program is free software: you can redistribute it and/or modify         //
// it under the terms of the version 3 of the GNU General Public License        //
// as published by the Free Software Foundation.                                //
//                                                                              //
// This program is distributed in the hope that it will be useful, but          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             //
// General Public License for more details.                                     //
//                                                                              //
// You should have received a copy of the GNU General Public License            //
// along with this program. If not, see <http://www.gnu.org/licenses/>.         //
//                                                                              //
// Contact <fethallah@gmail.com> for comments & bug reports                     //
//////////////////////////////////////////////////////////////////////////////////

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkSymmetricEigenAnalysis.h>
#include <itkImageRegionIterator.h>

const unsigned int maxDimension = 3;

#define SORTING_BY_ANISOTROPY_RATIO 1
#define SORTING_BY_MAGNITUDE        2


// Get the PixelType, the ComponentType and the number of dimensions
// from the fileName
void GetImageType (std::string							fileName,
                   itk::ImageIOBase::IOPixelType&		pixelType,
                   itk::ImageIOBase::IOComponentType&	componentType,
                   unsigned int&						noOfDimensions,
                   unsigned int&                        noOfComponents)
{
	typedef itk::Image<unsigned char, maxDimension>		ImageType;
	typedef itk::ImageFileReader<ImageType>				ReaderType;
	
	auto reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->UpdateOutputInformation();
	
	pixelType       = reader->GetImageIO()->GetPixelType();
	componentType   = reader->GetImageIO()->GetComponentType();
	noOfDimensions  = reader->GetImageIO()->GetNumberOfDimensions();
    noOfComponents  = reader->GetImageIO()->GetNumberOfComponents();
    
}

template<class TSymmetricDefinitePositveImageType >
void writeRiemannianMetricToTensorGlyphFile(typename TSymmetricDefinitePositveImageType::Pointer outputMetricImage, std::string glyphFileName, double percent)
{
    // TODO TODO TODO
    // first make the list of locations and the associated list of tensors
    // one could first sort all the pixels in outputMetricImage according to:
    //          1/ their anisotropy ratio
    //          2/ their maximal speed (which is the sqrt of the highest eigenvalue)
    //          3/ or just make a random selection
    // once they are sorted(selected), one can select a proportion of the best ones, for instance the 10% highest
    // TODO TODO TODO
        
    unsigned int const Dimension = TSymmetricDefinitePositveImageType::ImageDimension;
    typedef typename TSymmetricDefinitePositveImageType::PointType      LocationType;
    typedef typename TSymmetricDefinitePositveImageType::PixelType      TensorType;
    typedef itk::Vector<double, Dimension>                              EigenValuesType;
    typedef itk::Matrix<double, Dimension>                              EigenVectorsType;
    typedef itk::SymmetricEigenAnalysis
    <TensorType,
    EigenValuesType,
    EigenVectorsType>                                                   EigenAnalysisType;

    std::vector<LocationType>   listOfLocations;    listOfLocations.clear();
    std::vector<TensorType>     listOfTensors;      listOfTensors.clear();
    std::vector<double>         listOfScores;       listOfTensors.clear();
    
    EigenAnalysisType symmetricEigenAnalysis(Dimension);
    EigenValuesType  eigenValues;
    
    auto tensorsIt = itk::ImageRegionConstIterator<TSymmetricDefinitePositveImageType> (outputMetricImage, outputMetricImage->GetBufferedRegion());
    
    // first travers all the tensors, get the max speed along a given direction and sort them accordingly
    std::vector<double> listOfMaxSpeeds; listOfMaxSpeeds.clear();
    tensorsIt.GoToBegin();
    while (!tensorsIt.IsAtEnd())
    {
        symmetricEigenAnalysis.ComputeEigenValues(tensorsIt.Get(), eigenValues );
        listOfMaxSpeeds.push_back(eigenValues[Dimension-1]);
        ++tensorsIt;
    }
    std::sort(listOfMaxSpeeds.begin(), listOfMaxSpeeds.end()) ;
    
    unsigned int idxPercent = floor(listOfMaxSpeeds.size() *(1-percent));
    double threshold_value = listOfMaxSpeeds[idxPercent];
    
    tensorsIt.GoToBegin();
    while (!tensorsIt.IsAtEnd())
    {
        TensorType currentTensor = tensorsIt.Get();
        symmetricEigenAnalysis.ComputeEigenValues(currentTensor, eigenValues );
        LocationType pt;
        if(eigenValues[Dimension-1] > threshold_value)
        {
            auto idx = tensorsIt.GetIndex();
            outputMetricImage->TransformIndexToPhysicalPoint(idx, pt);
            listOfLocations.push_back(pt);
            listOfTensors.push_back(currentTensor);
            listOfScores.push_back(eigenValues[Dimension-1]);
        }
        ++tensorsIt;
    }
        
    // write to a file
    FILE* fid = fopen(glyphFileName.c_str(), "w");
	//ASCII file header
	fprintf(fid, "# vtk DataFile Version 3.0\n");
	fprintf(fid, "Random data to test tensors \n");
	fprintf(fid, "ASCII\n");
	fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fid, "POINTS %d float\n", int(listOfLocations.size()));
	fclose(fid);
	//append binary x,y,z data
	fid = fopen(glyphFileName.c_str(), "a");
	for (unsigned int i = 0; i < listOfLocations.size(); i++)
	{
		for (unsigned int j = 0; j < Dimension; j++)
		{
			float value = float( (listOfLocations[i])[j] );
			fprintf(fid, "%f ", value);
		}
        if(Dimension == 2)
        {
            fprintf(fid, "%f ", 0.0);
        }
		fprintf(fid, "\n");
	}
	    
    fprintf(fid, "POINT_DATA %d\n", int(listOfLocations.size()));
    fprintf(fid, "SCALARS scalars float\n");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for (unsigned int i = 0; i < listOfLocations.size(); i++)
	{
		float value = listOfScores[i];
		fprintf(fid, "%f ", value);
	}
    fprintf(fid, "\n");
        
	fprintf(fid, "\nTENSORS %s float\n", "RiemannianMetric");
	for (unsigned int i = 0; i < listOfLocations.size(); i++)
	{
        
        symmetricEigenAnalysis.ComputeEigenValues(listOfTensors[i], eigenValues );
		for (unsigned int j = 0; j < Dimension; j++)
		{
            for(unsigned int k = 0; k < Dimension; k++)
            {
                float value = (listOfTensors[i])(j, k) / eigenValues[Dimension-1]; // to normalize the size
                fprintf(fid, "%f ", value);
            }
            if(Dimension == 2)
            {
                fprintf(fid, "%f ", 0.0);
            }
            fprintf(fid, "\n");
            
		}
        if(Dimension == 2)
        {
            fprintf(fid, "0.0 0.0 0.0\n");
        }
        
		fprintf(fid, "\n");
	}
	fclose(fid);
}


template<class TInputImageType >
int GenerateSpatialRiemannianMetricFromOrientedFluxMatrix(int argc, char* argv[])
/**
 This case is very simple
 */
{
    const unsigned int                              Dimension = TInputImageType::ImageDimension;
    typedef typename TInputImageType::PixelType     SymmetrixMatrixType;
    typedef itk::Vector<double, Dimension>          EigenValuesType;
    typedef itk::Matrix<double, Dimension>          EigenVectorsType;
    
    typedef itk::SymmetricEigenAnalysis
        <SymmetrixMatrixType,
        EigenValuesType,
        EigenVectorsType>                           EigenAnalysisType;
    // Parse the input arguments.
	unsigned int argumentOffset = 1;
	std::string inputOFImageFilePath            = argv[argumentOffset++];
    std::string outputRiemannianMetricFilePath  = argv[argumentOffset++];
    bool        brightObjects                   = (bool)atoi(argv[argumentOffset++]);
    double      maxAnisotropyRatio              = atof(argv[argumentOffset++]);
    double      scaleSpeedRatio                 = atof(argv[argumentOffset++]);
    bool        generateGlyph                   = (bool)atoi(argv[argumentOffset++]);
    std::string glyphFileName                   = argv[argumentOffset++];
    double percentVisu                          = atof(argv[argumentOffset++]);
    
    auto OFMatrixImageReader = itk::ImageFileReader<TInputImageType>::New();
    OFMatrixImageReader->SetFileName(inputOFImageFilePath);
    try
    {
        OFMatrixImageReader->Update();
    }
    catch (itk::ExceptionObject &ex)
    {
        std::cout << ex << std::endl;
		return EXIT_FAILURE;
    }
    
    auto OFMatrixImage = OFMatrixImageReader->GetOutput();
//    OFMatrixImage->DisconnectPipeline(); Fucking buggy !!!
    
    auto outputMetricImage = TInputImageType::New();
    outputMetricImage->CopyInformation(OFMatrixImage);
    outputMetricImage->SetBufferedRegion(OFMatrixImage->GetLargestPossibleRegion());
    outputMetricImage->Allocate();
    
    EigenAnalysisType symmetricEigenAnalysis(Dimension);
    EigenValuesType  eigenValues;
    EigenVectorsType eigenVectors;
    // first get max_Lambda_n_minus_Lambda_1
    double max_Lambda_n_minus_Lambda_1 = 0;
    
    auto inIt  = itk::ImageRegionIterator<TInputImageType>(OFMatrixImage, OFMatrixImage->GetLargestPossibleRegion());
    inIt.GoToBegin();
    while(!inIt.IsAtEnd())
    {
        symmetricEigenAnalysis.ComputeEigenValues(inIt.Get(), eigenValues );
        if( max_Lambda_n_minus_Lambda_1 < (eigenValues[Dimension-1] - eigenValues[0]) )
        {
            max_Lambda_n_minus_Lambda_1 = (eigenValues[Dimension-1] - eigenValues[0]);
        }
        ++inIt;
    }
    
    itkAssertOrThrowMacro(max_Lambda_n_minus_Lambda_1 > 0, "max_Lambda_n_minus_Lambda_1 = 0, not interesting image of computation went wrong");
    // compute alpha
    double alpha;
    if(brightObjects)
    {
        alpha = 2*log(maxAnisotropyRatio) / max_Lambda_n_minus_Lambda_1;
    }
    else
    {
        alpha = -2*log(maxAnisotropyRatio) / max_Lambda_n_minus_Lambda_1;
    }
    
    auto outIt = itk::ImageRegionIterator<TInputImageType>(outputMetricImage, OFMatrixImage->GetLargestPossibleRegion());
    inIt.GoToBegin();
    outIt.GoToBegin();
    
    while(!inIt.IsAtEnd())
    {
        SymmetrixMatrixType outputMetric;
        outputMetric.Fill(0.0);
        symmetricEigenAnalysis.ComputeEigenValuesAndVectors(inIt.Get(), eigenValues, eigenVectors );
        
        for(unsigned int k = 0; k < Dimension; k++)
        {
            double sum_lambda_without_k = 0.0;
            for(unsigned int l = 0; l < Dimension; l++)
            {
                if(l != k)
                {
                    sum_lambda_without_k += eigenValues[l];
                }
            }
            
            for(unsigned int i = 0; i < Dimension; i++)
            {
                for(unsigned int j = 0; j < Dimension; j++)
                {
                    outputMetric(i, j) += exp(-alpha*sum_lambda_without_k)*eigenVectors[k][i]*eigenVectors[k][j];
                }
            }
        }
        
        outIt.Set(outputMetric);
        
        ++inIt;
        ++outIt;
    }
    
    auto outputMetricWriter = itk::ImageFileWriter<TInputImageType>::New();
    outputMetricWriter->SetFileName(outputRiemannianMetricFilePath);
    outputMetricWriter->SetInput(outputMetricImage);
    
    try
    {
        outputMetricWriter->Update();
    }
    catch (itk::ExceptionObject &ex)
    {
        std::cout << ex << std::endl;
		return EXIT_FAILURE;
    }
    
    if(generateGlyph)
    {
        writeRiemannianMetricToTensorGlyphFile<TInputImageType>(outputMetricImage, glyphFileName, percentVisu);
    }
    return EXIT_SUCCESS;
}


template<class TInputImageType>
int GenerateScaleSpaceRiemannianMetricFromOrientedFluxMatrix(int argc, char* argv[])
{
    return EXIT_SUCCESS;
}


void Usage(char* argv[])
{
	std::cerr << "Usage: " << std::endl
	<< argv[0] << std::endl
	<< " <input Oriented Flux Image Matrix Image file (could be scale space or space only)> " << std::endl
	<< " <output Riemannian Metric file> " << std::endl
    << " <BrightOjects> " << std::endl
    << " <maximum anisotropy ratio> " << std::endl
    << " <scale Speed Ratio (needed only for the scale space case)> " << std::endl
    << " <generate glyph (for visulizing the Riemannian metric)> " << std::endl
    << " <if the previous is true, where to write the file to be visualized with Paraview> " << std::endl
    << " <percentage of best tensors to visualize according to the sorting by magnitude> " << std::endl
	<< std::endl << std::endl;
}


// Check the arguments and try to parse the input image.
int main ( int argc, char* argv[] )
{
	if(argc < 9)
	{
		Usage(argv);
		return EXIT_FAILURE;
	}
	
	itk::ImageIOBase::IOPixelType		pixelType;
	itk::ImageIOBase::IOComponentType	componentType;
	unsigned int						noOfDimensions;
    unsigned int						noOfComponents;
	
	try
	{
		GetImageType(argv[1], pixelType, componentType, noOfDimensions, noOfComponents);
        
        if(pixelType == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR)
        {
            std::cout << noOfDimensions << std::endl;
            switch( noOfDimensions )
            {
                case 2:
                    switch (componentType)
                    {
                        case itk::ImageIOBase::DOUBLE:
                            typedef itk::SymmetricSecondRankTensor<double, 2>   SymMatrixSpaceDoubleType;
                            itkAssertOrThrowMacro(noOfComponents == 3, "nb components not appropriate");
                            GenerateSpatialRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixSpaceDoubleType, 2 > >(argc, argv);
                            break;
                        case itk::ImageIOBase::FLOAT:
                            typedef itk::SymmetricSecondRankTensor<float, 2>    SymMatrixSpaceFloatType;
                            itkAssertOrThrowMacro(noOfComponents == 3, "nb components not appropriate");
                            GenerateSpatialRiemannianMetricFromOrientedFluxMatrix< itk::Image<SymMatrixSpaceFloatType, 2> >(argc, argv);
                            break;
                        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                        default:
                            std::cout << "Unknown pixel component type: admissible types are float or double" << std::endl;
                            break;
                    }
                    break;
                
                case 3:
                    switch (componentType)
                    {
                        case itk::ImageIOBase::DOUBLE:
                            typedef itk::SymmetricSecondRankTensor<double, 3>   SymMatrixSpaceDoubleType;
                            typedef itk::SymmetricSecondRankTensor<double, 2>   SymMatrixScaleSpaceDoubleType;
                            if(noOfComponents == 3)
                            {
                                GenerateScaleSpaceRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixScaleSpaceDoubleType, 3> >(argc, argv);
                            }
                            else if(noOfComponents == 6)
                            {
                                GenerateSpatialRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixSpaceDoubleType, 3> >(argc, argv);
                            }
                            else
                            {
                                itkAssertOrThrowMacro(noOfComponents != 3 && noOfComponents != 6 , "nb components not appropriate");
                            }
                            break;
                        case itk::ImageIOBase::FLOAT:
                            typedef itk::SymmetricSecondRankTensor<float, 3>    SymMatrixSpaceFloatType;
                            typedef itk::SymmetricSecondRankTensor<float, 2>    SymMatrixScaleSpaceFloatType;
                            if(noOfComponents == 3)
                            {
                                GenerateScaleSpaceRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixScaleSpaceFloatType, 3> >(argc, argv);
                            }
                            else if(noOfComponents == 6)
                            {
                                GenerateSpatialRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixSpaceFloatType, 3> >(argc, argv);
                            }
                            else
                            {
                                itkAssertOrThrowMacro(noOfComponents != 3 && noOfComponents != 6 , "nb components not appropriate");
                            }
                            break;
                        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                        default:
                            std::cout << "Unknown pixel component type: admissible types are float or double" << std::endl;
                            break;
                    }
                    break;
                case 4:
                    switch (componentType)
                    {
                        case itk::ImageIOBase::DOUBLE:
                            typedef itk::SymmetricSecondRankTensor<double, 3>   SymMatrixScaleSpaceDoubleType;
                            itkAssertOrThrowMacro(noOfComponents == 6, "nb components not appropriate");
                            GenerateScaleSpaceRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixScaleSpaceDoubleType, 3> >(argc, argv);
                            break;
                        case itk::ImageIOBase::FLOAT:
                            typedef itk::SymmetricSecondRankTensor<float, 3>    SymMatrixScaleSpaceFloatType;
                            itkAssertOrThrowMacro(noOfComponents == 6, "nb components not appropriate");
                            GenerateScaleSpaceRiemannianMetricFromOrientedFluxMatrix<itk::Image<SymMatrixScaleSpaceFloatType, 3> >(argc, argv);
                            break;
                        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                        default:
                            std::cout << "Unknown pixel component type: admissible types are float or double" << std::endl;
                            break;

                    }
                    break;
                default:
                    std::cout << std::endl;
                    std::cout << "ERROR: possibile image dimensions are:" << std::endl;
                    std::cout << "       2: for space only metrics in 2D space" << std::endl;
                    std::cout << "       3: for space only metrics in 3D space, or scale space metrics in 2D space" << std::endl;
                    std::cout << "       4: for scale space metrics in 3D space" << std::endl;
                    break;
            }
        }
        else
        {
            std::cout << std::endl;
            std::cout << "ERROR: admissible input pixel type is SymmetricSecondRankTensor" << std::endl;
            return EXIT_FAILURE;
        }
	}
	catch( itk::ExceptionObject &excep)
	{
		std::cerr << argv[0] << ": exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


//symmetricEigenAnalysis.ComputeEigenValuesAndVectors(inIt.Get(), eigenValues, eigenVectors );

/**********************************************************/
/* Testin Code : tests passed successfully */
/**********************************************************/
/*
 itkAssertOrThrowMacro(eigenValues[0] <= eigenValues[1], "Eigen values not properly sorted");
 if(Dimension ==3)
 {
 itkAssertOrThrowMacro(eigenValues[1] <= eigenValues[2], "Eigen values not properly sorted");
 }
 for(unsigned int i = 0; i < Dimension; i++)
 {
 double normEigVector = 0;
 EigenValuesType EigenVector;
 for(unsigned int j = 0; j < Dimension; j++)
 {
 normEigVector += eigenVectors[i][j]*eigenVectors[i][j];
 EigenVector[j] = eigenVectors[i][j];
 }
 
 EigenValuesType ProdEigenVector;
 ProdEigenVector = 0.0 * ProdEigenVector;
 for(unsigned int j = 0; j < Dimension; j++)
 {
 for(unsigned int k = 0; k < Dimension; k++)
 {
 ProdEigenVector[j] += inIt.Get()(j, k)*EigenVector[k];
 }
 }
 ProdEigenVector = ProdEigenVector - eigenValues[i] * EigenVector;
 double normDiff = ProdEigenVector.GetNorm();
 itkAssertOrThrowMacro(normDiff < 1e-20, "norm of an eigen vector is not 1");
 itkAssertOrThrowMacro(fabs(normEigVector - 1) < 1e-6, "norm of an eigen vector is not 1");
 }
 */
