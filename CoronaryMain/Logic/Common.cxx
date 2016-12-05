#include "Common.h"
#include "LearningImpl.h"
#include "vtkSlicerCoronaryMainLogic.h"
#include <omp.h>



void FillIntegralImage(vtkImageData* intergalImage, vtkImageData *imageData, vtkImageInterpolator* interpolator)
{
	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];

	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);

	intergalImage->SetExtent(0, (int)(imageDims[0] * imageSpacings[0]) + 72, 0, (int)(imageDims[1] * imageSpacings[1]) + 72, 0, (int)(imageDims[2] * imageSpacings[2]) + 72);
	intergalImage->SetOrigin(imageOrigins[0] - 36.0, imageOrigins[1] - 36.0, imageOrigins[2] - 36.0);
	intergalImage->SetSpacing(1.0, 1.0, 1.0);
	intergalImage->AllocateScalars(VTK_DOUBLE, 6);

	int* dims = intergalImage->GetDimensions();
	double* origin = intergalImage->GetOrigin();
	double* spacing = intergalImage->GetSpacing();
	double coord[3];
	int dim[3];
	int pim[7][3];
	double pixel;
	int contri;
	int sign[][4] = { { -1, 0, 0, 1 }, { 0, -1, 0, 1 }, { -1, -1, 0, -1 }, { 0, 0, -1, 1 }, { -1, 0, -1, -1 }, { 0, -1, -1, -1 }, { -1, -1, -1, 1 } };

	double *integralimage = static_cast<double *>(intergalImage->GetScalarPointer());
	int dims01 = dims[0] * dims[1] * 6;
	int dims0 = dims[0] * 6;

	for (dim[2] = 0; dim[2] < dims[2]; dim[2]++)
	{
		for (dim[1] = 0; dim[1] < dims[1]; dim[1]++)
		{
			for (dim[0] = 0; dim[0] < dims[0]; dim[0]++)
			{
				for (int k = 0; k < 3; k++) coord[k] = origin[k] + spacing[k] * dim[k];
				interpolator->Interpolate(coord, &pixel);
				if (pixel < -512.0) contri = 0;
				else if (pixel < 0.0) contri = 1;
				else if (pixel < 256.0) contri = 2;
				else if (pixel < 512.0) contri = 3;
				else if (pixel < 768.0) contri = 4;
				else				  contri = 5;
				double* currhist = &integralimage[dim[2] * dims01 + dim[1] * dims0 + dim[0] * 6];
				memset(currhist, 0, 6 * sizeof(double));
				currhist[contri] = 1.0;

				for (int j = 0; j < 7; j++)
				{
					bool neg = false;
					for (int k = 0; k < 3; k++)
					{
						pim[j][k] = dim[k] + sign[j][k];
						neg = pim[j][k] < 0;
						if (neg) break;
					}
					if (neg) continue;

					double* prevhist = &integralimage[pim[j][2] * dims01 + pim[j][1] * dims0 + pim[j][0] * 6];
					for (int i = 0; i < 6; i++)
					{
						currhist[i] += sign[j][3] * prevhist[i];
					}
				}
			}
		}
	}
}

void IntegralImageHist(vtkImageData* intergalImage, int corner1[3], int corner2[3], double hist[6])
{
	int dim[3];
	double* inthist;
	for (int i = 0; i < 6; i++) hist[i] = 0.0;
	int sign[][4] = { { 0, 0, 0, 1 }, { -1, 0, 0, -1 }, { 0, -1, 0, -1 }, { -1, -1, 0, 1 }, { 0, 0, -1, -1 }, { -1, 0, -1, 1 }, { 0, -1, -1, 1 }, { -1, -1, -1, -1 } };
	for (int d = 0; d < 8; d++)
	{
		dim[0] = sign[d][0] == 0 ? corner2[0] : corner1[0];
		dim[1] = sign[d][1] == 0 ? corner2[1] : corner1[1];
		dim[2] = sign[d][2] == 0 ? corner2[2] : corner1[2];
		inthist = static_cast<double*>(intergalImage->GetScalarPointer(dim));
		for (int i = 0; i < 6; i++) hist[i] += sign[d][3] * inthist[i];
	}
}

void ImageFeatures(vtkImageData* intergalImage, double coord[3], cv::Mat& featureRow)
{
	double* origin = intergalImage->GetOrigin();
	double* spacing = intergalImage->GetSpacing();

	int dim[3];
	for (int i = 0; i < 3; i++) dim[i] = (int)((coord[i] - origin[i]) / spacing[i]);
	int corner1[3], corner2[3];

	double hist1[6], hist2[6];
	int count = 0;
	for (int i = 1; i <= 3; i++)
	{
		int cellsize = i * 8;

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2];
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize * 3 / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] - cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] + cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize * 3 / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize * 3 / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] - cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] + cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize * 3 / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);

		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize * 3 / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] - cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] - cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist2, hist1, hist1);
		corner1[0] = dim[0] - cellsize / 2; corner1[1] = dim[1] - cellsize / 2; corner1[2] = dim[2] + cellsize / 2;
		corner2[0] = dim[0] + cellsize / 2; corner2[1] = dim[1] + cellsize / 2; corner2[2] = dim[2] + cellsize * 3 / 2;
		IntegralImageHist(intergalImage, corner1, corner2, hist2);
		HistSubtract(hist1, hist2, hist1);
		InsertFeature(featureRow, hist1, count);
	}
}

void HistSubtract(const double hist1[6], const double hist2[6], double hist3[6])
{
	for (int i = 0; i < 6; i++) hist3[i] = hist1[i] - hist2[i];
}

void InsertFeature(cv::Mat& featureRow, double hist[6], int& count)
{
	HistNormalize(hist);
	for (int i = 0; i < 6; i++)
	{
		featureRow.at<float>(count) = hist[i];
		count++;
	}
}

void HistNormalize(double hist[6])
{
	double mag = 0.0;
	for (int i = 0; i < 6; i++) mag += hist[i] * hist[i];
	mag = std::max(sqrt(mag), 0.001);
	for (int i = 0; i < 6; i++) hist[i] /= mag;
}

void GenerateHessianImage(vtkImageData *imageData, vtkImageData *hessianImage, vtkImageInterpolator *interpolator)
{
	vtkSmartPointer<vtkImageResample> resample = vtkSmartPointer<vtkImageResample>::New();
	resample->SetInputData(imageData);
	for (int k = 0; k < 3; k++) resample->SetAxisMagnificationFactor(k, 0.5);
	resample->Update();
	vtkImageData* resampleImage = resample->GetOutput();

	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];
	resampleImage->GetDimensions(imageDims);
	resampleImage->GetOrigin(imageOrigins);
	resampleImage->GetSpacing(imageSpacings);

	vtkSmartPointer<vtkImageData> sumImage = vtkSmartPointer<vtkImageData>::New();
	double sumSpacing[3] = { 0.5, 0.5, 0.5 };
	double sumOrigin[3] = { imageOrigins[0] - 20.0, imageOrigins[1] - 20.0, imageOrigins[2] - 20.0 };
	int	   sumExtent[3];
	for (int k = 0; k < 3; k++) sumExtent[k] = int((imageDims[k] * imageSpacings[k] + 40.0) / sumSpacing[k]);

	sumImage->SetExtent(0, sumExtent[0], 0, sumExtent[1], 0, sumExtent[2]);
	sumImage->SetOrigin(sumOrigin);
	sumImage->SetSpacing(sumSpacing);
	sumImage->AllocateScalars(VTK_DOUBLE, 1);

	FillSumImage(sumImage, interpolator);

	//progressBar->setValue(40);

	hessianImage->CopyStructure(resampleImage);
	hessianImage->AllocateScalars(VTK_DOUBLE, 1);

	double coord[3];
	double eigvalue[3];
	double eigvector[3][3];

	short  *resampleimage = static_cast<short*>(resampleImage->GetScalarPointer());
	double *sumimage = static_cast<double*>(sumImage->GetScalarPointer());
	double *hessianimage = static_cast<double*>(hessianImage->GetScalarPointer());

	int dims01 = imageDims[1] * imageDims[0];
	int dims0 = imageDims[0];
	int scale[] = { 0, 7, 5, 3, 1 };
	for (int s = 0; s < 5; s++)
	{
#pragma omp parallel private(coord, eigvalue, eigvector) shared(dims01, dims0)
		{
#pragma omp for
			for (int d2 = 0; d2 < imageDims[2]; d2++)
			{
				int tid = omp_get_thread_num();
				if (tid == 0)
				{
					int myend = (tid + 1)*imageDims[2] / omp_get_num_threads();
				}

				for (int d1 = 0; d1 < imageDims[1]; d1++)
				{
					for (int d0 = 0; d0 < imageDims[0]; d0++)
					{
						if (s == 0)
						{
							hessianimage[d2*dims01 + d1*dims0 + d0] = 0.0;
							continue;
						}
						short pixel = resampleimage[d2*dims01 + d1*dims0 + d0];
						if (pixel < 100) continue;

						coord[0] = imageOrigins[0] + imageSpacings[0] * d0;
						coord[1] = imageOrigins[1] + imageSpacings[1] * d1;
						coord[2] = imageOrigins[2] + imageSpacings[2] * d2;

						Hessian(sumImage, sumimage, coord, eigvalue, eigvector, scale[s]);

						if (eigvalue[1] < 0.0)
						{
							//double RA2 = (eigvalue[1]*eigvalue[1])/(eigvalue[2]*eigvalue[2]);
							//double RB2 = (eigvalue[0]*eigvalue[0])/abs(eigvalue[1]*eigvalue[2]);
							//double RC2 = eigvalue[0]*eigvalue[0]+eigvalue[1]*eigvalue[1]+eigvalue[2]*eigvalue[2];
							//hessianimage[d2*dims01+d1*dims0+d0] = (1.0-exp(-RA2/0.5)) * exp(-RB2/0.5) * (1-exp(-RC2/1000.0));
							double lineMeasure;
							double normalizeValue = -eigvalue[1];
							if (eigvalue[0] < 0)
							{
								lineMeasure = eigvalue[0] / (0.5 * normalizeValue);
							}
							else
							{
								lineMeasure = eigvalue[0] / (2.0 * normalizeValue);
							}
							lineMeasure = exp(-0.5 * lineMeasure*lineMeasure);
							lineMeasure *= normalizeValue * eigvalue[1] / eigvalue[2];
							if (lineMeasure > hessianimage[d2*dims01 + d1*dims0 + d0])
								hessianimage[d2*dims01 + d1*dims0 + d0] = lineMeasure;
						}
						//else
						//{
						//	hessianimage[d2*dims01+d1*dims0+d0] = 0.0;
						//}
					}
				}
			}
		}
	}

	{
		typedef itk::Image<short, 3>  ImageType;
		typedef itk::VTKImageToImageFilter<ImageType> ResampleToImageType;
		ResampleToImageType::Pointer resampleToImageFilter = ResampleToImageType::New();
		resampleToImageFilter->SetInput(resampleImage);
		resampleToImageFilter->Update();

		typedef itk::Image<char, 3>  BinaryImageType;
		typedef itk::BinaryThresholdImageFilter<ImageType, BinaryImageType> thresholdFilterType;
		thresholdFilterType::Pointer thresholdFilter = thresholdFilterType::New();
		thresholdFilter->SetLowerThreshold(-2000);
		thresholdFilter->SetUpperThreshold(-250);
		thresholdFilter->SetInsideValue(1);
		thresholdFilter->SetOutsideValue(0);
		thresholdFilter->SetInput(resampleToImageFilter->GetOutput());

		typedef itk::BinaryBallStructuringElement<BinaryImageType::PixelType, BinaryImageType::ImageDimension> StructuringElementType;
		StructuringElementType structuringElement;
		structuringElement.SetRadius(5);
		structuringElement.CreateStructuringElement();

		typedef itk::BinaryMorphologicalClosingImageFilter<BinaryImageType, BinaryImageType, StructuringElementType> BinaryMorphologicalClosingImageFilterType;
		BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
		closingFilter->SetInput(thresholdFilter->GetOutput());
		closingFilter->SetKernel(structuringElement);
		closingFilter->SetForegroundValue(1);
		closingFilter->Update();
		//SaveITKImage<BinaryImageType>(closingFilter->GetOutput(), "lungseg.mha");

		typedef itk::ImageRegionConstIterator<BinaryImageType> BRegionIterator;
		BRegionIterator	rIt(closingFilter->GetOutput(), closingFilter->GetOutput()->GetRequestedRegion());
		BinaryImageType::IndexType index;
		for (rIt.GoToBegin(); !rIt.IsAtEnd(); ++rIt)
		{
			if (rIt.Get() > 0)
			{
				index = rIt.GetIndex();
				hessianimage[index[2] * dims01 + index[1] * dims0 + index[0]] = 0.0;
			}
		}
	}

	{
		typedef itk::Image<short, 3>  ImageType;
		typedef itk::VTKImageToImageFilter<ImageType> ResampleToImageType;
		ResampleToImageType::Pointer resampleToImageFilter = ResampleToImageType::New();
		resampleToImageFilter->SetInput(resampleImage);
		resampleToImageFilter->Update();

		typedef itk::Image<char, 3>  BinaryImageType;
		typedef itk::BinaryThresholdImageFilter<ImageType, BinaryImageType> thresholdFilterType;
		thresholdFilterType::Pointer thresholdFilter = thresholdFilterType::New();
		thresholdFilter->SetLowerThreshold(200);
		thresholdFilter->SetUpperThreshold(2000);
		thresholdFilter->SetInsideValue(1);
		thresholdFilter->SetOutsideValue(0);
		thresholdFilter->SetInput(resampleToImageFilter->GetOutput());

		typedef itk::BinaryBallStructuringElement<BinaryImageType::PixelType, BinaryImageType::ImageDimension> StructuringElementType;
		StructuringElementType structuringElement;
		structuringElement.SetRadius(8);
		structuringElement.CreateStructuringElement();

		typedef itk::BinaryMorphologicalOpeningImageFilter<BinaryImageType, BinaryImageType, StructuringElementType> BinaryMorphologicalOpeningImageFilterType;
		BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter = BinaryMorphologicalOpeningImageFilterType::New();
		openingFilter->SetInput(thresholdFilter->GetOutput());
		openingFilter->SetKernel(structuringElement);
		openingFilter->SetForegroundValue(1);
		openingFilter->Update();
		//SaveITKImage<BinaryImageType>(openingFilter->GetOutput(), "bloodpool.mha");

		typedef itk::ImageRegionConstIterator<BinaryImageType> BRegionIterator;
		BRegionIterator	rIt(openingFilter->GetOutput(), openingFilter->GetOutput()->GetRequestedRegion());
		BinaryImageType::IndexType index;
		for (rIt.GoToBegin(); !rIt.IsAtEnd(); ++rIt)
		{
			if (rIt.Get() > 0)
			{
				index = rIt.GetIndex();
				hessianimage[index[2] * dims01 + index[1] * dims0 + index[0]] = 0.0;
			}
		}
	}
}

void FillSumImage(vtkImageData* sumImage, vtkImageInterpolator* interpolator)
{
	int* dims = sumImage->GetDimensions();
	double* origin = sumImage->GetOrigin();
	double* spacing = sumImage->GetSpacing();
	double coord[3];
	int dim[3];
	int pim[3];
	double pixel;
	int sign[][4] = { { -1, 0, 0, 1 }, { 0, -1, 0, 1 }, { -1, -1, 0, -1 }, { 0, 0, -1, 1 }, { -1, 0, -1, -1 }, { 0, -1, -1, -1 }, { -1, -1, -1, 1 } };

	int dims01 = dims[0] * dims[1];
	double *image = static_cast<double*>(sumImage->GetScalarPointer());

	for (dim[2] = 0; dim[2] < dims[2]; dim[2]++)
	{
		for (dim[1] = 0; dim[1] < dims[1]; dim[1]++)
		{
			for (dim[0] = 0; dim[0] < dims[0]; dim[0]++)
			{
				for (int k = 0; k < 3; k++) coord[k] = origin[k] + spacing[k] * dim[k];
				interpolator->Interpolate(coord, &pixel);
				double& curr = image[dim[2] * dims01 + dim[1] * dims[0] + dim[0]];
				curr = pixel;

				for (int j = 0; j < 7; j++)
				{
					bool neg = false;
					for (int k = 0; k < 3; k++)
					{
						pim[k] = dim[k] + sign[j][k];
						neg = pim[k] < 0;
						if (neg) break;
					}
					if (neg) continue;

					curr += sign[j][3] * image[pim[2] * dims01 + pim[1] * dims[0] + pim[0]];
				}
			}
		}
	}
}

/* allocate memory for an nrow x ncol matrix */
template<class TReal>
TReal **create_matrix(long nrow, long ncol)
{
	typedef TReal* TRealPointer;
	TReal **m = new TRealPointer[nrow];

	TReal* block = (TReal*)calloc(nrow*ncol, sizeof(TReal));
	m[0] = block;
	for (int row = 1; row < nrow; ++row)
	{
		m[row] = &block[row * ncol];
	}
	return m;
}

/* free a TReal matrix allocated with matrix() */
template<class TReal>
void free_matrix(TReal **m)
{
	free(m[0]);
	delete[] m;
}


void SumImageHist(vtkImageData* sumImage, double* sumimage, int corner1[3], int corner2[3], double& sum)
{
	int* dims = sumImage->GetDimensions();
	int dims01 = dims[0] * dims[1];

	int dim[3];
	sum = 0.0;
	int sign[][4] = { { 0, 0, 0, 1 }, { -1, 0, 0, -1 }, { 0, -1, 0, -1 }, { -1, -1, 0, 1 }, { 0, 0, -1, -1 }, { -1, 0, -1, 1 }, { 0, -1, -1, 1 }, { -1, -1, -1, -1 } };
	for (int d = 0; d < 8; d++)
	{
		dim[0] = sign[d][0] == 0 ? corner2[0] : (corner1[0] - 1);
		dim[1] = sign[d][1] == 0 ? corner2[1] : (corner1[1] - 1);
		dim[2] = sign[d][2] == 0 ? corner2[2] : (corner1[2] - 1);
		sum += sign[d][3] * sumimage[dim[2] * dims01 + dim[1] * dims[0] + dim[0]];
	}
}



void Hessian(vtkImageData* sumImage, double* sumimage, double coord[3], double eigvalue[3], double eigvector[3][3], int cellsize = 3)
{
	double* origin = sumImage->GetOrigin();
	double* spacing = sumImage->GetSpacing();

	int dim[3];
	for (int i = 0; i < 3; i++) dim[i] = (int)((coord[i] - origin[i]) / spacing[i] + 0.5);
	int corner1[3], corner2[3];

	int cellsize3 = cellsize*cellsize*cellsize;
	int cellsizem3d2 = cellsize * 3 / 2;
	int cellsized2 = cellsize / 2;

	double sum1, sum2;
	double **hessian = create_matrix<double>(3, 3);
	//double hessian[3][3];
	{
		//dx^2
		corner1[0] = dim[0] - cellsizem3d2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] - cellsized2 - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] + cellsized2 + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsizem3d2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][0] = (sum1 + sum2) / cellsize3;

		//dy^2
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsizem3d2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] - cellsized2 - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] + cellsized2 + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsizem3d2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[1][1] = (sum1 + sum2) / cellsize3;

		//dz^2
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsizem3d2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] - cellsized2 - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - 2.0*sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] + cellsized2 + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsizem3d2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[2][2] = (sum1 + sum2) / cellsize3;

		//dxdy
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0]; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0]; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsized2;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsized2;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][1] = (sum1 - sum2) / cellsize3;
		hessian[1][0] = hessian[0][1];

		//dxdz
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsize + 1; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2];
		corner2[0] = dim[0]; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0]; corner1[1] = dim[1] - cellsized2; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsize - 1; corner2[1] = dim[1] + cellsized2; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[0][2] = (sum1 - sum2) / cellsize3;
		hessian[2][0] = hessian[0][2];

		//dydz
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1]; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum1);
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1]; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 + sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1] - cellsize + 1; corner1[2] = dim[2];
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1]; corner2[2] = dim[2] + cellsize - 1;
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		sum1 = sum1 - sum2;
		corner1[0] = dim[0] - cellsized2; corner1[1] = dim[1]; corner1[2] = dim[2] - cellsize + 1;
		corner2[0] = dim[0] + cellsized2; corner2[1] = dim[1] + cellsize - 1; corner2[2] = dim[2];
		SumImageHist(sumImage, sumimage, corner1, corner2, sum2);
		hessian[1][2] = (sum1 - sum2) / cellsize3;
		hessian[2][1] = hessian[1][2];
	}
	double **eigvec = create_matrix<double>(3, 3);
	vtkMath::Jacobi(hessian, eigvalue, eigvec);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			eigvector[i][j] = eigvec[i][j];

	free_matrix(hessian);
	free_matrix(eigvec);
}

template<class BinaryImageType>
void TraceCenterlineInternal(std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> >& clqueue,
	typename BinaryImageType::Pointer rawCenterline, vtkPoints *clPoints, vtkPolyData* clModel)
{
	std::map<typename BinaryImageType::IndexType, vtkIdType, typename BinaryImageType::IndexType::LexicographicCompare> indexMap;
	typedef itk::NeighborhoodIterator< BinaryImageType > NeighborhoodIteratorType;
	typename NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(1);
	NeighborhoodIteratorType nIt(nRadius, rawCenterline, rawCenterline->GetRequestedRegion());

	bool inbounds;

	while (!clqueue.empty())
	{
		std::list<typename BinaryImageType::IndexType> indexList;

		typename BinaryImageType::IndexType startid = clqueue.front().first;
		indexList.push_back(startid);
		typename BinaryImageType::IndexType nextid = clqueue.front().second;
		clqueue.pop();

		std::queue<typename BinaryImageType::IndexType> sgqueue;
		sgqueue.push(nextid);
		while (!sgqueue.empty())
		{
			nextid = sgqueue.front();
			sgqueue.pop();
			indexList.push_back(nextid);
			nIt.SetLocation(nextid);
			for (size_t i = 0; i < nIt.Size(); i++)
			{
				if (nIt.GetPixel(i, inbounds) > 0 && inbounds)
				{
					nIt.SetPixel(i, 0);
					sgqueue.push(nIt.GetIndex(i));
				}
			}
			if (sgqueue.size() != 1)
			{
				//Only insert when segment is long enough or segment is critical to the connectivity of the centerline
				if (indexList.size() > 10 || !sgqueue.empty())
				{
					//std::cout << "  cell " << clModel->GetNumberOfCells() << "  indexList: " << indexList.size() << std::endl;
					vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
					for (auto lit = indexList.begin(); lit != indexList.end(); ++lit)
					{
						if (indexMap.find(*lit) == indexMap.end())
						{
							typename BinaryImageType::PointType point;
							rawCenterline->TransformIndexToPhysicalPoint(*lit, point);
							indexMap[*lit] = clPoints->InsertNextPoint(point[0], point[1], point[2]);
						}
						idlist->InsertNextId(indexMap[*lit]);
					}
					clModel->InsertNextCell(VTK_POLY_LINE, idlist);
				}

				while (!sgqueue.empty())
				{
					clqueue.push(std::make_pair(nextid, sgqueue.front()));
					sgqueue.pop();
				}
			}
		}
	}
}

template<class BinaryImageType>
void TraceCenterline(typename BinaryImageType::Pointer rawCenterline, const typename BinaryImageType::IndexType& ostiumIndex, vtkPolyData* clModel)
{
	//std::cout << "Entering TraceCenterline: " << ostiumIndex << std::endl;

	vtkSmartPointer<vtkPoints> clPoints = vtkSmartPointer<vtkPoints>::New();

	typedef itk::NeighborhoodIterator< BinaryImageType > NeighborhoodIteratorType;
	NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(1);
	NeighborhoodIteratorType nIt(nRadius, rawCenterline, rawCenterline->GetRequestedRegion());

	bool inbounds;

	if (ostiumIndex[0] >= 0 && ostiumIndex[1] >= 0 && ostiumIndex[2] >= 0)		//find ostium seeded centerlines
	{
		nIt.SetLocation(ostiumIndex);
		if (nIt.GetCenterPixel() <= 0)	return;
		nIt.SetCenterPixel(-1);

		std::vector<typename BinaryImageType::IndexType> initdirs;
		std::vector<int>									 initsizes;
		std::vector< std::vector<typename BinaryImageType::IndexType> > initids;
		for (size_t i = 0; i < nIt.Size(); i++)
		{
			if (nIt.GetPixel(i, inbounds) > 0 && inbounds)
			{
				initdirs.push_back(nIt.GetIndex(i));
				initsizes.push_back(0);
				initids.push_back(std::vector<typename BinaryImageType::IndexType>());
			}
		}

		for (size_t i = 0; i < initdirs.size(); i++)
		{
			std::queue<typename BinaryImageType::IndexType> queue;
			queue.push(initdirs[i]);
			while (!queue.empty())
			{
				typename BinaryImageType::IndexType pid = queue.front();
				queue.pop();
				nIt.SetLocation(pid);

				nIt.SetCenterPixel(-1);
				initsizes[i]++;
				initids[i].push_back(pid);
				for (size_t i = 0; i < nIt.Size(); i++)
				{
					if (nIt.GetPixel(i, inbounds) > 0 && inbounds) queue.push(nIt.GetIndex(i));
				}
			}
		}

		auto maxdir = std::distance(initsizes.begin(), std::max_element(initsizes.begin(), initsizes.end()));
		for (size_t i = 0; i < initdirs.size(); i++)
		{
			if (i == maxdir)
			{
				for (size_t j = 0; j < initids[i].size(); j++)
				{
					nIt.SetLocation(initids[i][j]);
					nIt.SetCenterPixel(1);
				}
			}
			else
			{
				for (size_t j = 0; j < initids[i].size(); j++)
				{
					nIt.SetLocation(initids[i][j]);
					nIt.SetCenterPixel(0);
				}
			}
		}


		nIt.SetLocation(ostiumIndex);
		nIt.SetCenterPixel(0);
		nIt.SetLocation(initdirs[maxdir]);
		nIt.SetCenterPixel(0);
		std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> > clqueue;
		clqueue.push(std::make_pair(ostiumIndex, initdirs[maxdir]));


		TraceCenterlineInternal<BinaryImageType>(clqueue, rawCenterline, clPoints, clModel);
	}
	else
	{
		for (nIt.GoToBegin(); !nIt.IsAtEnd(); ++nIt)
		{
			if (nIt.GetCenterPixel() > 0)
			{
				for (size_t i = 0; i < nIt.Size(); i++)
				{
					if (nIt.GetPixel(i, inbounds) > 0 && inbounds)
					{
						nIt.SetLocation(nIt.GetIndex());
						nIt.SetCenterPixel(0);
						nIt.SetLocation(nIt.GetIndex(i));
						nIt.SetCenterPixel(0);
						std::queue< std::pair<typename BinaryImageType::IndexType, typename BinaryImageType::IndexType> > clqueue;
						clqueue.push(std::make_pair(nIt.GetIndex(), nIt.GetIndex(i)));
						TraceCenterlineInternal<BinaryImageType>(clqueue, rawCenterline, clPoints, clModel);
					}
				}
			}
		}
	}

	clModel->SetPoints(clPoints);
	//std::cout << "Leaving TraceCenterline: " << std::endl;
}

void CleanCenterline(vtkPolyData* clModel)
{
	if (!clModel || clModel->GetNumberOfCells() == 0) return;
	//std::cout << "Entering CleanCenterline:" << clModel->GetNumberOfCells() << std::endl;
	clModel->BuildCells();
	clModel->BuildLinks();

	unsigned short ncells;
	vtkIdType	*cells;
	vtkIdType numcells = clModel->GetNumberOfCells();
	std::vector<vtkIdType> deleteCells;
	for (vtkIdType i = 0; i < numcells; i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(i, idlist);
		vtkSmartPointer<vtkIdList> newlist = vtkSmartPointer<vtkIdList>::New();

		bool split = false;
		for (vtkIdType j = 0; j < idlist->GetNumberOfIds(); j++)
		{
			newlist->InsertNextId(idlist->GetId(j));

			clModel->GetPointCells(idlist->GetId(j), ncells, cells);
			if (ncells > 1 && j != 0 && j != idlist->GetNumberOfIds() - 1)
			{
				clModel->InsertNextCell(VTK_POLY_LINE, newlist);
				newlist->Initialize();
				newlist->InsertNextId(idlist->GetId(j));
				split = true;
			}
		}
		if (split)
		{
			clModel->InsertNextCell(VTK_POLY_LINE, newlist);
			deleteCells.push_back(i);
		}
	}
	if (deleteCells.size() > 0)
	{
		for (size_t i = 0; i < deleteCells.size(); i++)	clModel->DeleteCell(deleteCells[i]);
		clModel->RemoveDeletedCells();
		clModel->BuildCells();
		clModel->BuildLinks();
	}

	std::vector<int> newIds(clModel->GetNumberOfCells(), -1);
	int newid = 0;
	for (vtkIdType i = 0; i < clModel->GetNumberOfPoints(); i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetPointCells(i, idlist);
		if (idlist->GetNumberOfIds() == 2)
		{
			int& id1 = newIds[idlist->GetId(0)];
			int& id2 = newIds[idlist->GetId(1)];
			if (id1 < 0 && id2 < 0)
			{
				id1 = newid;
				id2 = newid;
				newid++;
			}
			else if (id1 >= 0 && id2 >= 0)
			{
				int oldid2 = id2, oldid1 = id1;
				if (oldid1 < oldid2)
				{

					for (int k = 0; k < newIds.size(); k++) if (newIds[k] == oldid2) newIds[k] = oldid1;
				}
				else if (oldid1 > oldid2)
				{
					for (int k = 0; k < newIds.size(); k++) if (newIds[k] == oldid1) newIds[k] = oldid2;
				}
			}
			else
			{
				if (id1 < id2) id1 = id2;
				else            id2 = id1;
			}
		}
	}
	vtkSmartPointer<vtkPolyData> newModel = vtkSmartPointer<vtkPolyData>::New();
	newModel->Allocate();
	for (vtkIdType i = 0; i < clModel->GetNumberOfCells(); i++)
	{
		if (newIds[i] == -1)
		{
			vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(i, idlist);
			newModel->InsertNextCell(VTK_POLY_LINE, idlist);
		}
		else if (newIds[i] >= 0)
		{
			int newid = newIds[i];
			std::set<vtkIdType> ids;

			for (int k = 0; k < newIds.size(); k++)
			{
				if (newIds[k] == newid)
				{
					ids.insert(k);
				}
			}
			//std::cout << "ids.size(): " << ids.size() << std::endl;

			std::deque<vtkIdType> pidlist;
			vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(*ids.begin(), idl);
			for (vtkIdType l = 0; l < idl->GetNumberOfIds(); l++) pidlist.push_back(idl->GetId(l));
			newIds[*ids.begin()] = -2;
			ids.erase(ids.begin());
			std::queue<vtkIdType> fqueue, equeue;
			fqueue.push(idl->GetId(0));
			equeue.push(idl->GetId(idl->GetNumberOfIds() - 1));
			while (!fqueue.empty())
			{
				vtkIdType fid = fqueue.front();
				fqueue.pop();
				for (auto iter = ids.begin(); iter != ids.end(); iter++)
				{
					vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(*iter, idl);
					if (idl->GetId(0) == fid)
					{
						for (vtkIdType l = 1; l < idl->GetNumberOfIds(); l++) pidlist.push_front(idl->GetId(l));
						fqueue.push(idl->GetId(idl->GetNumberOfIds() - 1));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
					else if (idl->GetId(idl->GetNumberOfIds() - 1) == fid)
					{
						for (vtkIdType l = idl->GetNumberOfIds() - 2; l >= 0; l--) pidlist.push_front(idl->GetId(l));
						fqueue.push(idl->GetId(0));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
				}
			}
			while (!equeue.empty())
			{
				vtkIdType fid = equeue.front();
				equeue.pop();
				for (auto iter = ids.begin(); iter != ids.end(); iter++)
				{
					vtkSmartPointer<vtkIdList> idl = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(*iter, idl);
					if (idl->GetId(0) == fid)
					{
						for (vtkIdType l = 1; l < idl->GetNumberOfIds(); l++) pidlist.push_back(idl->GetId(l));
						equeue.push(idl->GetId(idl->GetNumberOfIds() - 1));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
					else if (idl->GetId(idl->GetNumberOfIds() - 1) == fid)
					{
						for (vtkIdType l = idl->GetNumberOfIds() - 2; l >= 0; l--) pidlist.push_back(idl->GetId(l));
						equeue.push(idl->GetId(0));
						newIds[*iter] = -2;
						ids.erase(iter);
						break;
					}
				}
			}
			vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
			for (size_t p = 0; p < pidlist.size(); p++) idlist->InsertNextId(pidlist[p]);
			newModel->InsertNextCell(VTK_POLY_LINE, idlist);
		}
	}
	newModel->SetPoints(clModel->GetPoints());
	newModel->GetPointData()->CopyAllOn();
	newModel->GetPointData()->PassData(clModel->GetPointData());
	clModel->DeepCopy(newModel);

	//std::cout << "Leaving CleanCenterline:" << clModel->GetNumberOfCells() << std::endl;
}

void SimplifyCenterline(vtkPolyData* clModel)
{
	if (!clModel || clModel->GetNumberOfCells() == 0) return;
	//std::cout << "Entering SimplifyCenterline:" << clModel->GetNumberOfPoints() << std::endl;

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyData> newModel = vtkSmartPointer<vtkPolyData>::New();
	newModel->Allocate();
	std::map<vtkIdType, vtkIdType> pointMap;
	for (vtkIdType i = 0; i < clModel->GetNumberOfCells(); i++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(i, idlist);
		if (idlist->GetNumberOfIds() < 3) continue;

		std::list<double> polyline;
		double coord[3], last[3];
		double length = 0.0;
		for (vtkIdType l = 0; l < idlist->GetNumberOfIds(); l++)
		{
			clModel->GetPoint(idlist->GetId(l), coord);
			if (l > 0) length += sqrt(vtkMath::Distance2BetweenPoints(coord, last));
			for (int k = 0; k < 3; k++)
			{
				polyline.push_back(coord[k]);
				last[k] = coord[k];
			}
		}
		unsigned count = unsigned(length / 5.0);
		if (count<2)								count = 2;
		else if (count>idlist->GetNumberOfIds()) count = idlist->GetNumberOfIds();
		std::vector<double> newpolyline;
		psimpl::simplify_douglas_peucker_n<3>(polyline.begin(), polyline.end(), count, std::back_inserter(newpolyline));
		if (newpolyline.size() != 3 * count)
		{
			std::cerr << "Centerline simplification results wrong number of points" << std::endl;
			return;
		}

		double begin[3], end[3], newbegin[3], newend[3];
		clModel->GetPoint(idlist->GetId(0), begin);
		clModel->GetPoint(idlist->GetId(idlist->GetNumberOfIds() - 1), end);
		std::copy(newpolyline.begin(), newpolyline.begin() + 3, newbegin);
		std::copy(newpolyline.end() - 3, newpolyline.end(), newend);
		if (vtkMath::Distance2BetweenPoints(newbegin, begin) > 1e-6 || vtkMath::Distance2BetweenPoints(newend, end) > 1e-6)
		{
			std::cerr << "Centerline simplification results wrong begin and end points" << std::endl;
			return;
		}

		vtkSmartPointer<vtkIdList> newlist = vtkSmartPointer<vtkIdList>::New();
		for (size_t j = 0; j < count; j++)
		{
			for (int k = 0; k < 3; k++) coord[k] = newpolyline[3 * j + k];
			if (j == 0)
			{
				if (pointMap.find(idlist->GetId(0)) != pointMap.end())
				{
					newlist->InsertNextId(pointMap[idlist->GetId(0)]);
				}
				else
				{
					vtkIdType newid = newPoints->InsertNextPoint(coord);
					newlist->InsertNextId(newid);
					pointMap[idlist->GetId(0)] = newid;
				}
			}
			else if (j == count - 1)
			{
				if (pointMap.find(idlist->GetId(idlist->GetNumberOfIds() - 1)) != pointMap.end())
				{
					newlist->InsertNextId(pointMap[idlist->GetId(idlist->GetNumberOfIds() - 1)]);
				}
				else
				{
					vtkIdType newid = newPoints->InsertNextPoint(coord);
					newlist->InsertNextId(newid);
					pointMap[idlist->GetId(idlist->GetNumberOfIds() - 1)] = newid;
				}
			}
			else
			{
				vtkIdType newid = newPoints->InsertNextPoint(coord);
				newlist->InsertNextId(newid);
			}
		}
		newModel->InsertNextCell(VTK_POLY_LINE, newlist);
	}
	newModel->SetPoints(newPoints);

	vtkSmartPointer<ExtendSplineFilter>  clSpline = vtkSmartPointer<ExtendSplineFilter>::New();
	clSpline->SetSubdivideToLength();
	clSpline->SetLength(1.0);
	clSpline->SetGenerateTCoordsToOff();
	clSpline->SetInputData(newModel);
	clSpline->Update();

	clModel->DeepCopy(clSpline->GetOutput());

	//std::cout << "Leaving SimplifyCenterline:" << clModel->GetNumberOfPoints() << std::endl;
}

void RadiusCenterline(vtkPolyData* clModel)
{
	vtkSmartPointer<vtkDoubleArray> clRadius = vtkSmartPointer<vtkDoubleArray>::New();
	clRadius->SetName("Radius");
	clRadius->SetNumberOfValues(clModel->GetNumberOfPoints());
	for (vtkIdType id = 0; id < clModel->GetNumberOfPoints(); id++)
	{
		clRadius->SetValue(id, 1.0);
	}
	clModel->GetPointData()->SetScalars(clRadius);
}

void LumenWallCenterline(vtkPolyData* clModel)
{
	//std::cout << "Entering LumenWallCenterline" << std::endl;
	vtkDoubleArray* clRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Radius"));
	if (!clRadius) return;

	const int centerline_components = 12;	//must be even number
	vtkSmartPointer<vtkDoubleArray> clLumenRadius = vtkSmartPointer<vtkDoubleArray>::New();
	clLumenRadius->SetName("LumenRadius");
	clLumenRadius->SetNumberOfComponents(centerline_components);
	clLumenRadius->SetNumberOfTuples(clModel->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> clWallThickness = vtkSmartPointer<vtkDoubleArray>::New();
	clWallThickness->SetName("WallThickness");
	clWallThickness->SetNumberOfComponents(centerline_components);
	clWallThickness->SetNumberOfTuples(clModel->GetNumberOfPoints());

	double *radii = new double[clLumenRadius->GetNumberOfComponents()];
	for (vtkIdType id = 0; id < clModel->GetNumberOfPoints(); id++)
	{
		for (int j = 0; j < clLumenRadius->GetNumberOfComponents(); j++) radii[j] = clRadius->GetValue(id);
		clLumenRadius->SetTuple(id, radii);
		for (int j = 0; j < clWallThickness->GetNumberOfComponents(); j++) radii[j] = 0.2;
		clWallThickness->SetTuple(id, radii);
	}
	delete[] radii;
	clModel->GetPointData()->AddArray(clLumenRadius);
	clModel->GetPointData()->AddArray(clWallThickness);
	//std::cout << "Leaving LumenWallCenterline" << std::endl;
}


bool DetectLandmarks_core(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator)
{
	std::cout << "DetectLandmarks_core begin!" << std::endl;

	LearningImpl *learnimpl = learn.limpl;
	learnimpl->LoadLandmarkClassifiers(SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS);

	vtkSmartPointer<vtkImageData> integralImage = vtkSmartPointer<vtkImageData>::New();
	FillIntegralImage(integralImage, imageData, interpolator);
	std::cout << "FillIntegralImage done!" << std::endl;

	//SaveVTKImage(imageData, "C:\\work\\Coronary_Slicer\\testdata\\imageData.mha");
	//SaveVTKImage(integralImage, "C:\\work\\Coronary_Slicer\\testdata\\integralImage.mha");

	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];

	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);
	
	double coord[3];
	int dim[3];
	double maxpred[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS];
	for (int k = 0; k < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; k++) maxpred[k] = std::numeric_limits<double>::lowest();

	int imageDims01 = imageDims[0] * imageDims[1];
	short* imagedata = static_cast<short*>(imageData->GetScalarPointer());

	for (dim[2] = 10; dim[2] < imageDims[2] - 10; dim[2]++)
	{
		//if(dim[2]%10==0) progressBar->setValue(40+60*dim[2]/imageDims[2]);
		for (dim[1] = 20; dim[1] < imageDims[1] - 20; dim[1]++)
		{
			for (dim[0] = 20; dim[0] < imageDims[0] - 20; dim[0]++)
			{
				if ((dim[0] % 5 != 0 || dim[1] % 5 != 0 || dim[2] % 5 != 0)) continue;

				short pixel = imagedata[dim[2] * imageDims01 + dim[1] * imageDims[0] + dim[0]];
				if (pixel < 0) continue;

				for (int k = 0; k < 3; k++) coord[k] = imageOrigins[k] + imageSpacings[k] * dim[k];

				cv::Mat featureRow(1, 126, CV_32F);
				ImageFeatures(integralImage, coord, featureRow);
				for (int id = SmartCoronary::LEFT_CORONARY_OSTIUM; id < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; id++)
				{
					float pred = learnimpl->lmBoost[id].predict(featureRow, cv::Mat(), cv::Range::all(), false, true);
					if (pred > maxpred[id])
					{
						maxpred[id] = pred;
						for (int k = 0; k < 3; k++) landmarks[id][k] = coord[k];
					}
				}
			}
		}
	}

	std::cout << "DetectLandmarks_core done!" << std::endl;

	return true;
}

void GetRotationMatrix(double axis[3], double angle, double rot[3][3])
{
	double c = cos(angle);
	double s = sin(angle);
	rot[0][0] = c + axis[0] * axis[0] * (1.0 - c);
	rot[0][1] = axis[0] * axis[1] * (1.0 - c) - s*axis[2];
	rot[0][2] = axis[0] * axis[2] * (1.0 - c) + s*axis[1];
	rot[1][0] = axis[0] * axis[1] * (1.0 - c) + s*axis[2];
	rot[1][1] = c + axis[1] * axis[1] * (1.0 - c);
	rot[1][2] = axis[1] * axis[2] * (1.0 - c) - s*axis[0];
	rot[2][0] = axis[0] * axis[2] * (1.0 - c) - s*axis[1];
	rot[2][1] = axis[1] * axis[2] * (1.0 - c) + s*axis[0];
	rot[2][2] = c + axis[2] * axis[2] * (1.0 - c);
}

void AxisCenterline(vtkPolyData* clModel, double planenormal[3])
{
	//std::cout << "Entering AxisCenterline" << std::endl;
	if (planenormal) vtkMath::Normalize(planenormal);
	vtkSmartPointer<vtkDoubleArray>	clDir = vtkSmartPointer<vtkDoubleArray>::New();
	clDir->SetName("Dir");
	clDir->SetNumberOfComponents(3);
	clDir->SetNumberOfTuples(clModel->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray>	clAxis1 = vtkSmartPointer<vtkDoubleArray>::New();
	clAxis1->SetName("Axis1");
	clAxis1->SetNumberOfComponents(3);
	clAxis1->SetNumberOfTuples(clModel->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray>	clAxis2 = vtkSmartPointer<vtkDoubleArray>::New();
	clAxis2->SetName("Axis2");
	clAxis2->SetNumberOfComponents(3);
	clAxis2->SetNumberOfTuples(clModel->GetNumberOfPoints());
	double coord[3], dir[3], axis1[3], axis2[3], olddir[3], oldaxis1[3];
	double rot[3][3];
	for (vtkIdType id = 0; id < clModel->GetNumberOfCells(); id++)
	{
		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(id, idlist);
		for (vtkIdType j = 0; j < idlist->GetNumberOfIds(); j++)
		{
			if (j == 0)
			{
				clModel->GetPoint(idlist->GetId(j), coord);
				clModel->GetPoint(idlist->GetId(j + 1), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
				if (planenormal)
				{
					for (int k = 0; k < 3; k++) axis1[k] = planenormal[k];
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
				else
				{
					vtkMath::Perpendiculars(dir, axis1, axis2, 0.0);
				}
			}
			else if (j == idlist->GetNumberOfIds() - 1)
			{
				clModel->GetPoint(idlist->GetId(j - 1), coord);
				clModel->GetPoint(idlist->GetId(j), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
			}
			else
			{
				clModel->GetPoint(idlist->GetId(j - 1), coord);
				clModel->GetPoint(idlist->GetId(j + 1), dir);
				vtkMath::Subtract(dir, coord, dir);
				vtkMath::Normalize(dir);
			}
			if (j > 0)
			{
				vtkMath::Cross(olddir, dir, axis2);
				if (vtkMath::Norm(axis2) == 0.0)
				{
					for (int k = 0; k < 3; k++) axis1[k] = oldaxis1[k];
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
				else
				{
					vtkMath::Normalize(axis2);
					double angle = acos(vtkMath::Dot(olddir, dir));
					GetRotationMatrix(axis2, angle, rot);
					vtkMath::Multiply3x3(rot, oldaxis1, axis1);
					vtkMath::Normalize(axis1);
					vtkMath::Cross(dir, axis1, axis2);
					vtkMath::Normalize(axis2);
				}
			}
			clDir->SetTuple(idlist->GetId(j), dir);
			clAxis1->SetTuple(idlist->GetId(j), axis1);
			clAxis2->SetTuple(idlist->GetId(j), axis2);
			for (int k = 0; k < 3; k++)
			{
				olddir[k] = dir[k];
				oldaxis1[k] = axis1[k];

			}
		}
	}
	clModel->GetPointData()->AddArray(clDir);
	clModel->GetPointData()->AddArray(clAxis1);
	clModel->GetPointData()->AddArray(clAxis2);
	//std::cout << "Leaving AxisCenterline" << std::endl;
}

template<class ImageType>
void SaveITKImage(typename ImageType::Pointer image, const char* fileName)
{
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fileName);
	writer->SetInput(image);
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

bool DetectCenterline_core(vtkImageData *ImageData, vtkImageData *hessianImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3])
{
	typedef itk::Image<double, 3> ImageType;
	typedef itk::Image<short, 3>  BinaryImageType;
	typedef itk::BinaryThinningImageFilter3D<BinaryImageType, BinaryImageType> ThinningFilter;
	ThinningFilter::Pointer thinningFilter = ThinningFilter::New();
	
	std::cout << "leftOstium 1: " << leftOstium[0] << leftOstium[1] << leftOstium[2] << std::endl;
	std::cout << "rightOstium 1: " << rightOstium[0] << rightOstium[1] << rightOstium[2] << std::endl;
	ImageType::IndexType leftOstiumIndex, rightOstiumIndex;

	{
		typedef itk::VTKImageToImageFilter<ImageType> HessianToImageType;
		HessianToImageType::Pointer hessianToImageFilter = HessianToImageType::New();
		hessianToImageFilter->SetInput(hessianImage);
		hessianToImageFilter->Update();

		ImageType::Pointer vesselnessImage = hessianToImageFilter->GetOutput();

		ImageType::PointType leftOstiumPoint, rightOstiumPoint;
		for (int k = 0; k < 3; k++) leftOstiumPoint[k] = leftOstium[k];
		for (int k = 0; k < 3; k++) rightOstiumPoint[k] = rightOstium[k];

		vesselnessImage->TransformPhysicalPointToIndex(leftOstiumPoint, leftOstiumIndex);
		vesselnessImage->TransformPhysicalPointToIndex(rightOstiumPoint, rightOstiumIndex);
		std::cout << "leftOstiumIndex 1: " << leftOstiumIndex << std::endl;
		std::cout << "rightOstiumIndex 1: " << rightOstiumIndex << std::endl;
		typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
		NeighborhoodIteratorType::RadiusType nRadius; nRadius.Fill(4);
		NeighborhoodIteratorType nIt(nRadius, vesselnessImage, vesselnessImage->GetRequestedRegion());
		nIt.SetLocation(leftOstiumIndex);
		itk::OffsetValueType leftDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
		std::cout << "l1 nIt.Size() = " << nIt.Size() << std::endl;

		for (size_t i = 0; i < nIt.Size(); i++)
		{
			double pixel = nIt.GetPixel(i);
			if (pixel > 100.0)
			{
				NeighborhoodIteratorType::OffsetType offset = nIt.GetOffset(i);
				itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
				if (dist < leftDistanceVessness)
				{
					leftOstiumIndex = nIt.GetIndex(i);
					leftDistanceVessness = dist;
				}
			}
		}
		nIt.SetLocation(rightOstiumIndex);
		itk::OffsetValueType rightDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
		std::cout << "r1 nIt.Size() = " << nIt.Size() << std::endl;
		for (size_t i = 0; i < nIt.Size(); i++)
		{
			double pixel = nIt.GetPixel(i);
			if (pixel > 100.0)
			{
				NeighborhoodIteratorType::OffsetType offset = nIt.GetOffset(i);
				itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
				if (dist < rightDistanceVessness)
				{
					rightOstiumIndex = nIt.GetIndex(i);
					rightDistanceVessness = dist;
				}
			}
		}
		std::cout << "leftOstiumIndex 2: " << leftOstiumIndex << std::endl;
		std::cout << "rightOstiumIndex 2: " << rightOstiumIndex << std::endl;
		typedef itk::ConnectedThresholdImageFilter<ImageType, BinaryImageType> ThresholdFilterType;
		ThresholdFilterType::Pointer thresholdFilter1 = ThresholdFilterType::New();
		thresholdFilter1->SetInput(vesselnessImage);
		thresholdFilter1->SetConnectivity(ThresholdFilterType::FullConnectivity);
		thresholdFilter1->SetReplaceValue(1);
		if (leftDistanceVessness < std::numeric_limits<itk::OffsetValueType>::max() / 2)  thresholdFilter1->AddSeed(leftOstiumIndex);
		else
		{
			std::cerr << "Cannot find the starting point of the left coronary artery" << std::endl;
			return false;
		}
		if (rightDistanceVessness < std::numeric_limits<itk::OffsetValueType>::max() / 2) thresholdFilter1->AddSeed(rightOstiumIndex);
		else
		{
			std::cerr << "Cannot find the starting point of the right coronary artery" << std::endl;
			return false;
		}
		thresholdFilter1->SetUpper(1500.0); // 1500
		thresholdFilter1->SetLower(100.0);  // 100
		thresholdFilter1->Update();

	/*	typedef itk::CastImageFilter< BinaryImageType, ImageType > CastFilterType;
		CastFilterType::Pointer castFilter = CastFilterType::New();
		castFilter->SetInput(thresholdFilter1->GetOutput());
		castFilter->Update();
		typedef itk::ImageToVTKImageFilter<ImageType> vtkFrangiResultType;
		vtkFrangiResultType::Pointer vtkFrangiResult = vtkFrangiResultType::New();
		vtkFrangiResult->SetInput(castFilter->GetOutput());
		vtkFrangiResult->Update();
	*/
	//	SaveITKImage<BinaryImageType>(thresholdFilter1->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\thresholdFilter.mha");
	//	SaveVTKImage(vtkFrangiResult->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\vtkFrangiResult.mha");

		thinningFilter->SetInput(thresholdFilter1->GetOutput());
		thinningFilter->Update();

	//	SaveITKImage<BinaryImageType>(thinningFilter->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\thinningFilter.mha");

		typedef itk::NeighborhoodIterator<BinaryImageType>	BNeighborhoodIteratorType;
		nRadius.Fill(5);
		nIt.SetRadius(nRadius);
		BNeighborhoodIteratorType tIt(nRadius, thinningFilter->GetOutput(), thinningFilter->GetOutput()->GetRequestedRegion());
		tIt.SetLocation(leftOstiumIndex);
		nIt.SetLocation(leftOstiumIndex);
		leftDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
		std::cout << "l2 nIt.Size() = " << nIt.Size() << std::endl;
		std::cout << "l2 tIt.Size() = " << tIt.Size() << std::endl;
		for (size_t i = 0; i<nIt.Size(); i++)
		{
			double pixel = nIt.GetPixel(i);
			if (tIt.GetPixel(i) && pixel > 100.0)
			{
				NeighborhoodIteratorType::OffsetType offset = nIt.GetOffset(i);
				itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
				if (dist < leftDistanceVessness)
				{
					leftOstiumIndex = nIt.GetIndex(i);
					leftDistanceVessness = dist;
				}
			}
		}
		tIt.SetLocation(rightOstiumIndex);
		nIt.SetLocation(rightOstiumIndex);
		rightDistanceVessness = std::numeric_limits<itk::OffsetValueType>::max();
		std::cout << "r2 nIt.Size() = " << nIt.Size() << std::endl;
		std::cout << "r2 tIt.Size() = " << tIt.Size() << std::endl;
		for (size_t i = 0; i<nIt.Size(); i++)
		{
			double pixel = nIt.GetPixel(i);
			if (tIt.GetPixel(i) && pixel > 100.0)
			{
				NeighborhoodIteratorType::OffsetType offset = nIt.GetOffset(i);
				itk::OffsetValueType dist = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
				if (dist < rightDistanceVessness)
				{
					rightOstiumIndex = nIt.GetIndex(i);
					rightDistanceVessness = dist;
				}
			}
		}
		if (leftDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2)
		{
			std::cerr << "Cannot find the starting point of the left coronary artery after thinning" << std::endl;
			return false;
		}
		if (rightDistanceVessness > std::numeric_limits<itk::OffsetValueType>::max() / 2)
		{
			std::cerr << "Cannot find the starting point of the right coronary artery after thinning" << std::endl;
			return false;
		}
		std::cout << "leftOstiumIndex 3: " << leftOstiumIndex << std::endl;
		std::cout << "rightOstiumIndex 3: " << rightOstiumIndex << std::endl;

	}

	typedef itk::ImageDuplicator< BinaryImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(thinningFilter->GetOutput());
	duplicator->Update();
	BinaryImageType::Pointer rawCenterline = duplicator->GetModifiableOutput();
	
//	SaveITKImage<BinaryImageType>(rawCenterline, "C:\\work\\Coronary_Slicer\\testdata\\rawCenterline.mha");
	
	//ImageType::IndexType invalidIndex; invalidIndex.Fill(-1);
	//ImageType::IndexType ostiumIndex[3] = {leftOstiumIndex, rightOstiumIndex, invalidIndex};
	ImageType::IndexType ostiumIndex[2] = { leftOstiumIndex, rightOstiumIndex };
	std::cout << "ostiumIndex[2]: " << ostiumIndex[0] << ", " << ostiumIndex[1] << std::endl;
	vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
	for (int i = 1; i >= 0; i--)
	{
		vtkSmartPointer<vtkPolyData> clModel = vtkSmartPointer<vtkPolyData>::New();
		clModel->Allocate();
		TraceCenterline<BinaryImageType>(rawCenterline, ostiumIndex[i], clModel);
		CleanCenterline(clModel);
		SimplifyCenterline(clModel);
		RadiusCenterline(clModel);
		LumenWallCenterline(clModel);
		AxisCenterline(clModel);
		append->AddInputData(clModel);
	}

	append->Update();
	centerlineModel->DeepCopy(append->GetOutput());
	
	std::cout << "DetectCenterline_core done!" << std::endl;

	return true;
}



void SaveVTKImage(vtkImageData *image, const char* fileName)
{
	vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
	writer->SetFileName(fileName);
	writer->SetInputData(image);
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

void SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}