#include "Common.h"
#include "LearningImpl.h"
#include "vtkSlicerCoronaryMainLogic.h"
#include <omp.h>



CBifurcation::CBifurcation()
{
	this->CenterPID = 0;
	this->VesselID.resize(0);
	this->EndfacePID.resize(0);
	this->VesselDir.resize(0);
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


int smoothvtkpolydata(vtkPolyData* Poly, int iternum, int TYPE)
{
	vtkPoints* Points = Poly->GetPoints();
	vtkCellArray* Strips = Poly->GetPolys();
	if (Strips->GetNumberOfCells() == 0)
	{
		Strips = Poly->GetStrips();
	}

//	if (TYPE == 1)
//		Strips = Poly->GetPolys();
//	else
//		Strips = Poly->GetStrips();

	if (Strips->GetNumberOfCells() == 0)
	{
		std::cerr << "Cannot find cell data" << std::endl;
		return 0;
	}

	vtkIdType CellId = 0;
	vtkIdType npts = 0, *pts = NULL;

	int* adjcent = (int*)malloc(Points->GetNumberOfPoints() * 10 * sizeof(int));
	int* num_adjcent = (int*)malloc(Points->GetNumberOfPoints() * sizeof(int));

	for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid++)
	{
		num_adjcent[pid] = 0;
	}
	for (CellId = 0, Strips->InitTraversal(); Strips->GetNextCell(npts, pts); CellId++)
	{
		if (npts != 3)
		{
			std::cout << "not triangle, smooth cannot work!" << std::endl;
			return 0;
		}
		for (int i = 0; i < npts; i++)
		{
			int p[2];
			int pidx = 0;
			for (int k = 0; k < npts; k++)
			{
				if (k != i)
				{
					p[pidx] = k;
					pidx++;
				}
			}
			for (int l = 0; l < 2; l++)
			{
				bool find_pl_in_adjofptsi = false;
				for (int k = 0; k < num_adjcent[pts[i]]; k++)
				{
					if (adjcent[pts[i] * 10 + k] == pts[p[l]])
					{
						find_pl_in_adjofptsi = true;
						break;
					}
				}

				if (find_pl_in_adjofptsi == false)
				{
					adjcent[pts[i] * 10 + num_adjcent[pts[i]]] = pts[p[l]];
					num_adjcent[pts[i]] ++;
				}
			}
		}
	}


	// the smooth algorithm
	{
		vtkSmartPointer<vtkPoints> Points_orig = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> Points_last = vtkSmartPointer<vtkPoints>::New();
		Points_orig->DeepCopy(Points);

		double* b = (double*)malloc(Points->GetNumberOfPoints() * 3 * sizeof(double));
		double pi[3], qi[3], oi[3];
		const double alpha = 0.1, beta = 0.2;

		for (int iter = 0; iter < iternum; iter++)
		{
			Points_last->DeepCopy(Points);
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid++)
			{
				if (num_adjcent[pid] > 0)
				{
					pi[0] = 0;
					pi[1] = 0;
					pi[2] = 0;
					for (int j = 0; j < num_adjcent[pid]; j++)
					{
						Points_last->GetPoint(adjcent[pid * 10 + j], qi);
						pi[0] += qi[0];
						pi[1] += qi[1];
						pi[2] += qi[2];
					}
					pi[0] = pi[0] / num_adjcent[pid];
					pi[1] = pi[1] / num_adjcent[pid];
					pi[2] = pi[2] / num_adjcent[pid];
					Points->SetPoint(pid, pi);
					Points_orig->GetPoint(pid, oi);
					Points_last->GetPoint(pid, qi);

					b[pid * 3 + 0] = pi[0] - (alpha * oi[0] + (1.0 - alpha) * qi[0]);
					b[pid * 3 + 1] = pi[1] - (alpha * oi[1] + (1.0 - alpha) * qi[1]);
					b[pid * 3 + 2] = pi[2] - (alpha * oi[2] + (1.0 - alpha) * qi[2]);

				}
			}
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid++)
			{
				if (num_adjcent[pid] > 0)
				{
					double sumbj[3];
					sumbj[0] = 0;
					sumbj[1] = 0;
					sumbj[2] = 0;

					for (int j = 0; j < num_adjcent[pid]; j++)
					{
						sumbj[0] += b[adjcent[pid * 10 + j] * 3 + 0];
						sumbj[1] += b[adjcent[pid * 10 + j] * 3 + 1];
						sumbj[2] += b[adjcent[pid * 10 + j] * 3 + 2];
					}
					sumbj[0] = sumbj[0] / num_adjcent[pid];
					sumbj[1] = sumbj[1] / num_adjcent[pid];
					sumbj[2] = sumbj[2] / num_adjcent[pid];

					Points->GetPoint(pid, pi);

					pi[0] = pi[0] - (beta * b[pid * 3 + 0] + (1.0 - beta) * sumbj[0]);
					pi[1] = pi[1] - (beta * b[pid * 3 + 1] + (1.0 - beta) * sumbj[1]);
					pi[2] = pi[2] - (beta * b[pid * 3 + 2] + (1.0 - beta) * sumbj[2]);

					Points->SetPoint(pid, pi);
				}
			}
		}
		Points_orig->Reset();
		Points_last->Reset();
		free(b);
	}

	Poly->SetPoints(Points);
	return 1;
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

void GenerateHessianImage(vtkImageData *imageData, vtkImageData *hessianImage, vtkImageInterpolator *interpolator, QProgressBar* progressbar)
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

	//-------//
	progressbar->setValue(10);

	FillSumImage(sumImage, interpolator);

	//-------//
	progressbar->setValue(15);

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

	//-------//
	progressbar->setValue(30);

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

	//-------//
	progressbar->setValue(40);
	
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

bool DetectLandmarks_core(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator, QProgressBar* progressbar)
{
	std::cout << "DetectLandmarks_core begin!" << std::endl;

	LearningImpl *learnimpl = learn.limpl;
	learnimpl->LoadLandmarkClassifiers(SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS);

	//-------// 
	progressbar->setValue(10);

	vtkSmartPointer<vtkImageData> integralImage = vtkSmartPointer<vtkImageData>::New();
	FillIntegralImage(integralImage, imageData, interpolator);

	//-------// 
	progressbar->setValue(40);

	std::cout << "FillIntegralImage done!" << std::endl;

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
		//-------// 
		if (dim[2] % 10 == 0)	progressbar->setValue(40 + 50 * dim[2] / (imageDims[2] - 20));

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

bool DetectCenterline_core(vtkImageData *ImageData, vtkImageData *hessianImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3], QProgressBar* progressbar)
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

		//-------// 
		progressbar->setValue(60);

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
//		std::cout << "l1 nIt.Size() = " << nIt.Size() << std::endl;

		//-------// 
		progressbar->setValue(70);

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
//		std::cout << "r1 nIt.Size() = " << nIt.Size() << std::endl;

		//-------// 
		progressbar->setValue(80);

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

		//-------// 
		progressbar->setValue(85);

		typedef itk::ConnectedThresholdImageFilter<ImageType, BinaryImageType> ThresholdFilterType;
		ThresholdFilterType::Pointer thresholdFilter1 = ThresholdFilterType::New();
		thresholdFilter1->SetInput(vesselnessImage);
		thresholdFilter1->SetConnectivity(ThresholdFilterType::FullConnectivity);
		thresholdFilter1->SetReplaceValue(1);
		if (leftDistanceVessness < std::numeric_limits<itk::OffsetValueType>::max() / 2)  thresholdFilter1->AddSeed(leftOstiumIndex);
		else
		{
			std::cerr << "Cannot find the starting point of the left coronary artery" << std::endl;
			//-------// 
			progressbar->setValue(0);
			return false;
		}
		if (rightDistanceVessness < std::numeric_limits<itk::OffsetValueType>::max() / 2) thresholdFilter1->AddSeed(rightOstiumIndex);
		else
		{
			std::cerr << "Cannot find the starting point of the right coronary artery" << std::endl;
			//-------// 
			progressbar->setValue(0);
			return false;
		}
		thresholdFilter1->SetUpper(1500.0); // 1500
		thresholdFilter1->SetLower(100.0);  // 100
		thresholdFilter1->Update();

		//-------// 
		progressbar->setValue(90);


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
	//	std::cout << "l2 nIt.Size() = " << nIt.Size() << std::endl;
	//	std::cout << "l2 tIt.Size() = " << tIt.Size() << std::endl;
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
	//	std::cout << "r2 nIt.Size() = " << nIt.Size() << std::endl;
	//	std::cout << "r2 tIt.Size() = " << tIt.Size() << std::endl;
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

		//-------// 
		progressbar->setValue(95);

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


bool DetectCenterlineLumenWall_core(vtkPolyData* clModel, vtkIdType selectId, vtkImageInterpolator* interpolator, Learning &learn)
{
	LearningImpl *learnimpl = learn.limpl;
	
	//	learnimpl->LoadLumenWallClassifiers();
	learnimpl->LoadLumenAlongNormalsClassifiers();

	vtkSmartPointer<vtkDoubleArray> clAxis1 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis1"));
	vtkSmartPointer<vtkDoubleArray> clAxis2 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis2"));
	vtkSmartPointer<vtkDoubleArray> clLumenRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
	vtkSmartPointer<vtkDoubleArray> clWallThickness = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("WallThickness"));
	if (!clAxis1 || !clAxis2 || !clLumenRadius || !clWallThickness) return false;

	vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
	clModel->GetCellPoints(selectId, idlist);

	double center[3], coord[3][3], axis1[3], axis2[3], ray[3][3];
	double cirstep = 2.0*M_PI / clLumenRadius->GetNumberOfComponents();
	double* Radius = new double[clLumenRadius->GetNumberOfComponents()];
	double* Radius_temp = new double[clLumenRadius->GetNumberOfComponents()];
	double* Thickness = new double[clLumenRadius->GetNumberOfComponents()];
	double* Preds = new double[clLumenRadius->GetNumberOfComponents()];
	for (vtkIdType i = 0; i < idlist->GetNumberOfIds(); i++)
	{
		vtkIdType pid = idlist->GetId(i);
		clModel->GetPoint(pid, center);
		clAxis1->GetTuple(pid, axis1);
		clAxis2->GetTuple(pid, axis2);

		for (int k = 0; k < clLumenRadius->GetNumberOfComponents(); k++)
		{
			for (int l = 0; l < 3; l++) ray[0][l] = cos((k - 1)*cirstep)*axis1[l] + sin((k - 1)*cirstep)*axis2[l];
			for (int l = 0; l < 3; l++) ray[1][l] = cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l];
			for (int l = 0; l < 3; l++) ray[2][l] = cos((k + 1)*cirstep)*axis1[l] + sin((k + 1)*cirstep)*axis2[l];
			double maxpred = std::numeric_limits<double>::lowest();
			double maxradius, maxthickness;
			for (double radius = 0.8; radius < 2.5; radius = radius + 0.05)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int l = 0; l < 3; l++) coord[j][l] = center[l] + radius*ray[j][l];
				}

				double pred = 0;
				for (double thickness = 0.05; thickness < 0.25; thickness += 0.05)
				{
					for (int j = 0; j < 3; j++)
					{
						cv::Mat featureRow = cv::Mat::zeros(1, 25, CV_32F);
						RayFeatures(interpolator, coord[j], ray[j], thickness, featureRow);
						pred += learnimpl->lwBoost.predict(featureRow, cv::Mat(), cv::Range::all(), false, true);
					}
				}
				if (pred > maxpred)
				{
					maxpred = pred;
					maxradius = radius;
					//maxthickness = thickness;
					maxthickness = 0.05;
				}
			}
			Radius[k] = maxradius;
			Thickness[k] = maxthickness;
			Preds[k] = maxpred;
		}

		double center_new[3];
		clModel->GetPoint(pid, center);
		CentralizedThisContour(center, axis1, axis2, clLumenRadius->GetNumberOfComponents(), center_new, Radius, Thickness);
		clModel->GetPoints()->SetPoint(pid, center_new);

		for (int k = 0; k < clLumenRadius->GetNumberOfComponents(); k++)
		{
			clLumenRadius->SetComponent(pid, k, Radius[k]);
			clWallThickness->SetComponent(pid, k, Thickness[k]);
		}
	}

	delete[] Radius_temp;
	delete[] Radius;
	delete[] Thickness;
	delete[] Preds;
	return true;
}





int MergeAlgorithm(vector<CEndFace> endfaces, double bifurcationcenter[3], vector<CBifurcationTriangle>& triangles)
{
	int number_endfaces = endfaces.size();
	int number_ringpoints = endfaces[0].rx.size();

	vector<CBifurcationTriangle> triangles_orignal;

	for (int fid = 0; fid < number_endfaces; fid ++)
	{
		for (int pid = 0; pid < number_ringpoints; pid++)
		{
			int pid2 = (pid == number_ringpoints - 1) ? 0 : pid + 1;
			{
				CBifurcationTriangle thistrianglemesh;
				findConvexPoint(fid, pid, pid2, endfaces, &thistrianglemesh);
			
				triangles_orignal.push_back(thistrianglemesh);
			}
		}
	}


	int trianglesize_orignal = triangles_orignal.size();
	bool findthisedge = false;
	for (int tid = 0; tid < trianglesize_orignal; tid ++)
	{
		for (int pid = 0; pid < triangles_orignal[tid].EndFacePoint.size(); pid++)
		{
			int pid2 = (pid == triangles_orignal[tid].EndFacePoint.size() - 1) ? 0 : pid + 1;
			findthisedge = false;

			if (triangles_orignal[tid].EndFacePoint[pid].index[0] == triangles_orignal[tid].EndFacePoint[pid2].index[0])
				continue;

			for (int tid_inside = 0; tid_inside < triangles_orignal.size(); tid_inside++)
			{
				if (tid_inside == tid)	continue;

				for (int pid_inside = 0; pid_inside < triangles_orignal[tid_inside].EndFacePoint.size(); pid_inside++)
				{
					int pid2_inside = (pid_inside == triangles_orignal[tid_inside].EndFacePoint.size() - 1) ? 0 : pid_inside + 1;

					if (triangles_orignal[tid_inside].EndFacePoint[pid_inside].index[0] == triangles_orignal[tid_inside].EndFacePoint[pid2_inside].index[0])
						continue;

					if ((triangles_orignal[tid_inside].EndFacePoint[pid_inside].index[0] == triangles_orignal[tid].EndFacePoint[pid].index[0] && triangles_orignal[tid_inside].EndFacePoint[pid_inside].index[1] == triangles_orignal[tid].EndFacePoint[pid].index[1]
						&& triangles_orignal[tid_inside].EndFacePoint[pid2_inside].index[0] == triangles_orignal[tid].EndFacePoint[pid2].index[0] && triangles_orignal[tid_inside].EndFacePoint[pid2_inside].index[1] == triangles_orignal[tid].EndFacePoint[pid2].index[1])
						|| (triangles_orignal[tid_inside].EndFacePoint[pid_inside].index[0] == triangles_orignal[tid].EndFacePoint[pid2].index[0] && triangles_orignal[tid_inside].EndFacePoint[pid_inside].index[1] == triangles_orignal[tid].EndFacePoint[pid2].index[1]
						&& triangles_orignal[tid_inside].EndFacePoint[pid2_inside].index[0] == triangles_orignal[tid].EndFacePoint[pid].index[0] && triangles_orignal[tid_inside].EndFacePoint[pid2_inside].index[1] == triangles_orignal[tid].EndFacePoint[pid].index[1]))
					{
						findthisedge = true;
						break;
					}
				}
				if (findthisedge) break;
			}
			if (findthisedge) continue;

			int pid3 = 0;
			// notice 			for (i3 = 0; i3 < localPointNuminEachContour[k] - 1; i3 ++)
			for (pid3 = 0; pid3 < triangles_orignal[tid].EndFacePoint.size(); pid3++)
			{
				if (pid3 != pid && pid3 != pid2)
					break;
			}

			CBifurcationTriangle thistrianglemesh;

			findConvexPoint_fill_big_triangle_hole(triangles_orignal[tid].EndFacePoint[pid].index[0], triangles_orignal[tid].EndFacePoint[pid2].index[0], triangles_orignal[tid].EndFacePoint[pid3].index[0]
												 , triangles_orignal[tid].EndFacePoint[pid].index[1], triangles_orignal[tid].EndFacePoint[pid2].index[1], triangles_orignal[tid].EndFacePoint[pid3].index[1]
												 , endfaces, &thistrianglemesh);
	
			triangles_orignal.push_back(thistrianglemesh);
		}
	}
	

	// manually sub divide convex hull triangles
	for (int i = 0; i < triangles_orignal.size(); i++)
	{
		// triangles_no_subdivision[i].EndFacePoint.index[0];
		int j1, j2, j3;
		
		if (   triangles_orignal[i].EndFacePoint[0].index[0]
			== triangles_orignal[i].EndFacePoint[1].index[0])
		{
			j1 = 2;
			j2 = 0;
			j3 = 1;
		}
		else if (triangles_orignal[i].EndFacePoint[0].index[0]
			  == triangles_orignal[i].EndFacePoint[2].index[0])
		{
			j1 = 1;
			j2 = 0;
			j3 = 2;
		}
		else if (triangles_orignal[i].EndFacePoint[1].index[0]
			  == triangles_orignal[i].EndFacePoint[2].index[0])
		{
			j1 = 0;
			j2 = 1;
			j3 = 2;
		}
		else
		{
			j1 = -1;
			j2 = -1;
			j3 = -1;
		}

		if (j1 != -1)
		{
			double mid1[3], mid2[3];			// divide 1/2
			double sub1[3], sub2[3], sub3[3], sub4[3]; // divide 1/4

			for (int l = 0; l < 3; l++)
			{
				mid1[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j1].realcoord[l] + triangles_orignal[i].EndFacePoint[j2].realcoord[l]);
				mid2[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j1].realcoord[l] + triangles_orignal[i].EndFacePoint[j3].realcoord[l]);
			}
			for (int l = 0; l < 3; l++)
			{
				sub1[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j1].realcoord[l] + mid1[l]);
				sub2[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j1].realcoord[l] + mid2[l]);
				sub3[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j2].realcoord[l] + mid1[l]);
				sub4[l] = 0.5 * (triangles_orignal[i].EndFacePoint[j3].realcoord[l] + mid2[l]);
			}

			CBifurcationTriangle thistrianglemesh;

			// [i][j1] sub1 sub2
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[j1].index[0], triangles_orignal[i].EndFacePoint[j1].index[1], triangles_orignal[i].EndFacePoint[j1].coord[0], triangles_orignal[i].EndFacePoint[j1].coord[1], triangles_orignal[i].EndFacePoint[j1].coord[2], triangles_orignal[i].EndFacePoint[j1].realcoord[0], triangles_orignal[i].EndFacePoint[j1].realcoord[1], triangles_orignal[i].EndFacePoint[j1].realcoord[2]
				, -1, -1, sub1[0], sub1[1], sub1[2], sub1[0], sub1[1], sub1[2]
				, -1, -1, sub2[0], sub2[1], sub2[2], sub2[0], sub2[1], sub2[2]);
			triangles.push_back(thistrianglemesh);

			// midcoord1 sub1 sub2
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, sub1[0], sub1[1], sub1[2], sub1[0], sub1[1], sub1[2]
				, -1, -1, sub2[0], sub2[1], sub2[2], sub2[0], sub2[1], sub2[2]);
			triangles.push_back(thistrianglemesh);

			// midcoord1 midcoord2 sub2
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub2[0], sub2[1], sub2[2], sub2[0], sub2[1], sub2[2]);
			triangles.push_back(thistrianglemesh);

			// midcoord1 midcoord2 sub3
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub3[0], sub3[1], sub3[2], sub3[0], sub3[1], sub3[2]);
			triangles.push_back(thistrianglemesh);

			// midcoord2 sub3 sub4
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub3[0], sub3[1], sub3[2], sub3[0], sub3[1], sub3[2]
				, -1, -1, sub4[0], sub4[1], sub4[2], sub4[0], sub4[1], sub4[2]);
			triangles.push_back(thistrianglemesh);

			// j2 sub3 sub4
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[j2].index[0], triangles_orignal[i].EndFacePoint[j2].index[1], triangles_orignal[i].EndFacePoint[j2].coord[0], triangles_orignal[i].EndFacePoint[j2].coord[1], triangles_orignal[i].EndFacePoint[j2].coord[2], triangles_orignal[i].EndFacePoint[j2].realcoord[0], triangles_orignal[i].EndFacePoint[j2].realcoord[1], triangles_orignal[i].EndFacePoint[j2].realcoord[2]
				, -1, -1, sub3[0], sub3[1], sub3[2], sub3[0], sub3[1], sub3[2]
				, -1, -1, sub4[0], sub4[1], sub4[2], sub4[0], sub4[1], sub4[2]);
			triangles.push_back(thistrianglemesh);

			// j2 j3 sub4
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[j2].index[0], triangles_orignal[i].EndFacePoint[j2].index[1], triangles_orignal[i].EndFacePoint[j2].coord[0], triangles_orignal[i].EndFacePoint[j2].coord[1], triangles_orignal[i].EndFacePoint[j2].coord[2], triangles_orignal[i].EndFacePoint[j2].realcoord[0], triangles_orignal[i].EndFacePoint[j2].realcoord[1], triangles_orignal[i].EndFacePoint[j2].realcoord[2]
				, triangles_orignal[i].EndFacePoint[j3].index[0], triangles_orignal[i].EndFacePoint[j3].index[1], triangles_orignal[i].EndFacePoint[j3].coord[0], triangles_orignal[i].EndFacePoint[j3].coord[1], triangles_orignal[i].EndFacePoint[j3].coord[2], triangles_orignal[i].EndFacePoint[j3].realcoord[0], triangles_orignal[i].EndFacePoint[j3].realcoord[1], triangles_orignal[i].EndFacePoint[j3].realcoord[2]
				, -1, -1, sub4[0], sub4[1], sub4[2], sub4[0], sub4[1], sub4[2]);
			triangles.push_back(thistrianglemesh);
		}
		else
		{
			double mid1[3], mid2[3], mid3[3];
			double sub11[3], sub12[3], sub13[3];
			double sub21[3], sub22[3], sub23[3];
			double sub31[3], sub32[3], sub33[3];

			for (int l = 0; l < 3; l++)
			{
				mid1[l] = 0.5 * (triangles_orignal[i].EndFacePoint[0].realcoord[l] + triangles_orignal[i].EndFacePoint[1].realcoord[l]);
				mid2[l] = 0.5 * (triangles_orignal[i].EndFacePoint[1].realcoord[l] + triangles_orignal[i].EndFacePoint[2].realcoord[l]);
				mid3[l] = 0.5 * (triangles_orignal[i].EndFacePoint[2].realcoord[l] + triangles_orignal[i].EndFacePoint[0].realcoord[l]);
			}
			for (int l = 0; l < 3; l++)
			{
				sub11[l] = 0.5 * (triangles_orignal[i].EndFacePoint[0].realcoord[l] + mid1[l]);
				sub12[l] = 0.5 * (triangles_orignal[i].EndFacePoint[0].realcoord[l] + mid3[l]);
				sub13[l] = 0.5 * (mid1[l] + mid3[l]);
				sub21[l] = 0.5 * (triangles_orignal[i].EndFacePoint[1].realcoord[l] + mid1[l]);
				sub22[l] = 0.5 * (triangles_orignal[i].EndFacePoint[1].realcoord[l] + mid2[l]);
				sub23[l] = 0.5 * (mid1[l] + mid2[l]);
				sub31[l] = 0.5 * (triangles_orignal[i].EndFacePoint[2].realcoord[l] + mid2[l]);
				sub32[l] = 0.5 * (triangles_orignal[i].EndFacePoint[2].realcoord[l] + mid3[l]);
				sub33[l] = 0.5 * (mid2[l] + mid3[l]);
			}

			CBifurcationTriangle thistrianglemesh;

			// 0 sub11 sub12
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[0].index[0], triangles_orignal[i].EndFacePoint[0].index[1], triangles_orignal[i].EndFacePoint[0].coord[0], triangles_orignal[i].EndFacePoint[0].coord[1], triangles_orignal[i].EndFacePoint[0].coord[2], triangles_orignal[i].EndFacePoint[0].realcoord[0], triangles_orignal[i].EndFacePoint[0].realcoord[1], triangles_orignal[i].EndFacePoint[0].realcoord[2]
				, -1, -1, sub11[0], sub11[1], sub11[2], sub11[0], sub11[1], sub11[2]
				, -1, -1, sub12[0], sub12[1], sub12[2], sub12[0], sub12[1], sub12[2]);
			triangles.push_back(thistrianglemesh);

			// mid1 sub11 sub13
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, sub11[0], sub11[1], sub11[2], sub11[0], sub11[1], sub11[2]
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]);
			triangles.push_back(thistrianglemesh);

			// mid3 sub12 sub13
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid3[0], mid3[1], mid3[2], mid3[0], mid3[1], mid3[2]
				, -1, -1, sub12[0], sub12[1], sub12[2], sub12[0], sub12[1], sub12[2]
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]);
			triangles.push_back(thistrianglemesh);

			// sub11 sub12 sub13
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, sub11[0], sub11[1], sub11[2], sub11[0], sub11[1], sub11[2]
				, -1, -1, sub12[0], sub12[1], sub12[2], sub12[0], sub12[1], sub12[2]
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]);
			triangles.push_back(thistrianglemesh);

			// 1 sub21 sub22
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[1].index[0], triangles_orignal[i].EndFacePoint[1].index[1], triangles_orignal[i].EndFacePoint[1].coord[0], triangles_orignal[i].EndFacePoint[1].coord[1], triangles_orignal[i].EndFacePoint[1].coord[2], triangles_orignal[i].EndFacePoint[1].realcoord[0], triangles_orignal[i].EndFacePoint[1].realcoord[1], triangles_orignal[i].EndFacePoint[1].realcoord[2]
				, -1, -1, sub21[0], sub21[1], sub21[2], sub21[0], sub21[1], sub21[2]
				, -1, -1, sub22[0], sub22[1], sub22[2], sub22[0], sub22[1], sub22[2]);
			triangles.push_back(thistrianglemesh);

			// mid1 sub21 sub23
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, sub21[0], sub21[1], sub21[2], sub21[0], sub21[1], sub21[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]);
			triangles.push_back(thistrianglemesh);

			// mid2 sub22 sub23
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub22[0], sub22[1], sub22[2], sub22[0], sub22[1], sub22[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]);
			triangles.push_back(thistrianglemesh);

			// sub21 sub22 sub23
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, sub21[0], sub21[1], sub21[2], sub21[0], sub21[1], sub21[2]
				, -1, -1, sub22[0], sub22[1], sub22[2], sub22[0], sub22[1], sub22[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]);
			triangles.push_back(thistrianglemesh);

			// 2 sub31 sub32
			fillBifurcationTriangle(&thistrianglemesh
				, triangles_orignal[i].EndFacePoint[2].index[0], triangles_orignal[i].EndFacePoint[2].index[1], triangles_orignal[i].EndFacePoint[2].coord[0], triangles_orignal[i].EndFacePoint[2].coord[1], triangles_orignal[i].EndFacePoint[2].coord[2], triangles_orignal[i].EndFacePoint[2].realcoord[0], triangles_orignal[i].EndFacePoint[2].realcoord[1], triangles_orignal[i].EndFacePoint[2].realcoord[2]
				, -1, -1, sub31[0], sub31[1], sub31[2], sub31[0], sub31[1], sub31[2]
				, -1, -1, sub32[0], sub32[1], sub32[2], sub32[0], sub32[1], sub32[2]);
			triangles.push_back(thistrianglemesh);

			// mid2 sub31 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub31[0], sub31[1], sub31[2], sub31[0], sub31[1], sub31[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

			// mid3 sub32 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid3[0], mid3[1], mid3[2], mid3[0], mid3[1], mid3[2]
				, -1, -1, sub32[0], sub32[1], sub32[2], sub32[0], sub32[1], sub32[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

			// sub31 sub32 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, sub31[0], sub31[1], sub31[2], sub31[0], sub31[1], sub31[2]
				, -1, -1, sub32[0], sub32[1], sub32[2], sub32[0], sub32[1], sub32[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

			// mid1 sub23 sub13
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid1[0], mid1[1], mid1[2], mid1[0], mid1[1], mid1[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]);
			triangles.push_back(thistrianglemesh);

			// mid2 sub23 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid2[0], mid2[1], mid2[2], mid2[0], mid2[1], mid2[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

			// mid3 sub13 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, mid3[0], mid3[1], mid3[2], mid3[0], mid3[1], mid3[2]
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

			//sub13 sub23 sub33
			fillBifurcationTriangle(&thistrianglemesh
				, -1, -1, sub13[0], sub13[1], sub13[2], sub13[0], sub13[1], sub13[2]
				, -1, -1, sub23[0], sub23[1], sub23[2], sub23[0], sub23[1], sub23[2]
				, -1, -1, sub33[0], sub33[1], sub33[2], sub33[0], sub33[1], sub33[2]);
			triangles.push_back(thistrianglemesh);

		}
	}

//	std::cout << "triangles.size() = " << triangles.size() << endl;

//	std::cout << "==========================" << endl;
	
	return 0;
}

int findConvexPoint(int fid
	, int pid, int pid2
	, vector<CEndFace> endfaces
	, CBifurcationTriangle* trianglemesh_out)
{
	int NUM_SEGMENT = endfaces.size();
	int RING_SIZE = endfaces[0].rx.size();

	double x = endfaces[fid].x;
	double y = endfaces[fid].y;
	double z = endfaces[fid].z;

	double ax = endfaces[fid].rx[pid2] - endfaces[fid].rx[pid];
	double ay = endfaces[fid].ry[pid2] - endfaces[fid].ry[pid];
	double az = endfaces[fid].rz[pid2] - endfaces[fid].rz[pid];

	double temp = sqrt(ax * ax + ay * ay + az * az);
	ax = ax / temp;
	ay = ay / temp;
	az = az / temp;

	double p[3], q[3], qp[3];
	p[0] = endfaces[fid].rx[0] - x;
	p[1] = endfaces[fid].ry[0] - y;
	p[2] = endfaces[fid].rz[0] - z;
	q[0] = endfaces[fid].rx[4] - x;
	q[1] = endfaces[fid].ry[4] - y;
	q[2] = endfaces[fid].rz[4] - z;
	vtkMath::Cross(p, q, qp);
	double dx = qp[0];
	double dy = qp[1];
	double dz = qp[2];

	temp = sqrt(dx * dx + dy * dy + dz * dz);
	dx = dx / temp;
	dy = dy / temp;
	dz = dz / temp;

	double ox = endfaces[fid].rx[pid];
	double oy = endfaces[fid].ry[pid];
	double oz = endfaces[fid].rz[pid];

	double cx = endfaces[fid].rx[pid] - x;
	double cy = endfaces[fid].ry[pid] - y;
	double cz = endfaces[fid].rz[pid] - z;
	double adx = ay * dz - az * dy;
	double ady = az * dx - dz * ax;
	double adz = ax * dy - ay * dx;

	if (adx * cx + ady * cy + adz * cz < 0)
	{
		adx = -adx;
		ady = -ady;
		adz = -adz;
	}

	double maxTh = std::numeric_limits<double>::lowest();
	double th = 0.0;

	int IDl = 0;		// one convex point of line (j,j2) is 
	int IDj = 0;
	double outx, outy, outz;

	double px, py, pz;

	for (int l = 0; l < NUM_SEGMENT; l++)
	{
		if (l == fid)
			continue;

		for (int j = 0; j < RING_SIZE; j++)
		{
			px = endfaces[l].rx[j] - ox - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * ax;
			py = endfaces[l].ry[j] - oy - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * ay;
			pz = endfaces[l].rz[j] - oz - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * az;

			temp = sqrt(px * px + py * py + pz * pz);
			px = px / temp;
			py = py / temp;
			pz = pz / temp;

			th = px * adx + py * ady + pz * adz;

			if (th > maxTh)
			{
				maxTh = th;

				IDl = l;
				IDj = j;
				outx = endfaces[l].rx[j];
				outy = endfaces[l].ry[j];
				outz = endfaces[l].rz[j];
			}
		}
	}

	fillBifurcationTriangle(trianglemesh_out
		, fid, pid, endfaces[fid].rx[pid], endfaces[fid].ry[pid], endfaces[fid].rz[pid], endfaces[fid].realrx[pid], endfaces[fid].realry[pid], endfaces[fid].realrz[pid]
		, fid, pid2, endfaces[fid].rx[pid2], endfaces[fid].ry[pid2], endfaces[fid].rz[pid2], endfaces[fid].realrx[pid2], endfaces[fid].realry[pid2], endfaces[fid].realrz[pid2]
		, IDl, IDj, endfaces[IDl].rx[IDj], endfaces[IDl].ry[IDj], endfaces[IDl].rz[IDj], endfaces[IDl].realrx[IDj], endfaces[IDl].realry[IDj], endfaces[IDl].realrz[IDj]);


	return 0;
}

int findConvexPoint_fill_big_triangle_hole(int fid1, int fid2, int fid3
	, int pid1, int pid2, int pid3
	, vector<CEndFace> endfaces
	, CBifurcationTriangle* trianglemesh_out)
{
	int NUM_SEGMENT = endfaces.size();
	int RING_SIZE = endfaces[0].rx.size();

	double ax = endfaces[fid2].rx[pid2] - endfaces[fid1].rx[pid1];
	double ay = endfaces[fid2].ry[pid2] - endfaces[fid1].ry[pid1];
	double az = endfaces[fid2].rz[pid2] - endfaces[fid1].rz[pid1];
	
	double temp = sqrt(ax * ax + ay * ay + az * az);
	ax = ax / temp;
	ay = ay / temp;
	az = az / temp;

	double x = endfaces[fid3].rx[pid3];
	double y = endfaces[fid3].ry[pid3];
	double z = endfaces[fid3].rz[pid3];

	double p[3], q[3], qp[3];
	p[0] = endfaces[fid1].rx[0] - x;
	p[1] = endfaces[fid1].ry[0] - y;
	p[2] = endfaces[fid1].rz[0] - z;
	q[0] = endfaces[fid1].rx[4] - x;
	q[1] = endfaces[fid1].ry[4] - y;
	q[2] = endfaces[fid1].rz[4] - z;
	vtkMath::Cross(p, q, qp);
	double dx = qp[0];
	double dy = qp[1];
	double dz = qp[2];

	temp = sqrt(dx * dx + dy * dy + dz * dz);
	dx = dx / temp;
	dy = dy / temp;
	dz = dz / temp;

	double ox = endfaces[fid1].rx[pid1];
	double oy = endfaces[fid1].ry[pid1];
	double oz = endfaces[fid1].rz[pid1];

	double cx = endfaces[fid1].rx[pid1] - x;
	double cy = endfaces[fid1].ry[pid1] - y;
	double cz = endfaces[fid1].rz[pid1] - z;
	double adx = ay * dz - az * dy;
	double ady = az * dx - dz * ax;
	double adz = ax * dy - ay * dx;

	if (adx * cx + ady * cy + adz * cz < 0)
	{
		adx = -adx;
		ady = -ady;
		adz = -adz;
	}

	double maxTh = std::numeric_limits<double>::lowest();
	double th = 0.0;

	int IDl = 0;		// one convex point of line (j,j2) is 
	int IDj = 0;
	double outx, outy, outz;

	double px, py, pz;

	for (int l = 0; l < NUM_SEGMENT; l++)
	{
		if (l == fid1 || l == fid2)
			continue;

		for (int j = 0; j < RING_SIZE; j++)
		{
			px = endfaces[l].rx[j] - ox - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * ax;
			py = endfaces[l].ry[j] - oy - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * ay;
			pz = endfaces[l].rz[j] - oz - (ax * (endfaces[l].rx[j] - ox) + ay * (endfaces[l].ry[j] - oy) + az * (endfaces[l].rz[j] - oz)) * az;

			temp = sqrt(px * px + py * py + pz * pz);
			px = px / temp;
			py = py / temp;
			pz = pz / temp;

			th = px * adx + py * ady + pz * adz;

			if (th > maxTh)
			{
				maxTh = th;

				IDl = l;
				IDj = j;
				outx = endfaces[l].rx[j];
				outy = endfaces[l].ry[j];
				outz = endfaces[l].rz[j];
			}
		}
	}
	fillBifurcationTriangle(trianglemesh_out
		, fid1, pid1, endfaces[fid1].rx[pid1], endfaces[fid1].ry[pid1], endfaces[fid1].rz[pid1], endfaces[fid1].realrx[pid1], endfaces[fid1].realry[pid1], endfaces[fid1].realrz[pid1]
		, fid2, pid2, endfaces[fid2].rx[pid2], endfaces[fid2].ry[pid2], endfaces[fid2].rz[pid2], endfaces[fid2].realrx[pid2], endfaces[fid2].realry[pid2], endfaces[fid2].realrz[pid2]
		, IDl, IDj, endfaces[IDl].rx[IDj], endfaces[IDl].ry[IDj], endfaces[IDl].rz[IDj], endfaces[IDl].realrx[IDj], endfaces[IDl].realry[IDj], endfaces[IDl].realrz[IDj]);

	return 0;
}

void fillBifurcationTriangle(CBifurcationTriangle* t
	, vtkIdType p1_idx0, vtkIdType p1_idx1, double p1_rx, double p1_ry, double p1_rz, double p1_realrx, double p1_realry, double p1_realrz
	, vtkIdType p2_idx0, vtkIdType p2_idx1, double p2_rx, double p2_ry, double p2_rz, double p2_realrx, double p2_realry, double p2_realrz
	, vtkIdType p3_idx0, vtkIdType p3_idx1, double p3_rx, double p3_ry, double p3_rz, double p3_realrx, double p3_realry, double p3_realrz)
{
	t->EndFacePoint.clear();
	CEndFacePoint EndFacePoint;

	EndFacePoint.index[0] = p1_idx0;
	EndFacePoint.index[1] = p1_idx1;
	EndFacePoint.coord[0] = p1_rx;
	EndFacePoint.coord[1] = p1_ry;
	EndFacePoint.coord[2] = p1_rz;
	EndFacePoint.realcoord[0] = p1_realrx;
	EndFacePoint.realcoord[1] = p1_realry;
	EndFacePoint.realcoord[2] = p1_realrz;
	t->EndFacePoint.push_back(EndFacePoint);

	EndFacePoint.index[0] = p2_idx0;
	EndFacePoint.index[1] = p2_idx1;
	EndFacePoint.coord[0] = p2_rx;
	EndFacePoint.coord[1] = p2_ry;
	EndFacePoint.coord[2] = p2_rz;
	EndFacePoint.realcoord[0] = p2_realrx;
	EndFacePoint.realcoord[1] = p2_realry;
	EndFacePoint.realcoord[2] = p2_realrz;
	t->EndFacePoint.push_back(EndFacePoint);

	EndFacePoint.index[0] = p3_idx0;
	EndFacePoint.index[1] = p3_idx1;
	EndFacePoint.coord[0] = p3_rx;
	EndFacePoint.coord[1] = p3_ry;
	EndFacePoint.coord[2] = p3_rz;
	EndFacePoint.realcoord[0] = p3_realrx;
	EndFacePoint.realcoord[1] = p3_realry;
	EndFacePoint.realcoord[2] = p3_realrz;
	t->EndFacePoint.push_back(EndFacePoint);
}

void ImageGradient(vtkImageInterpolator* interpolator, const double point[3], double gradient[3], double scale = 1.0)
{
	double coord[3];
	double grad1, grad2;
	for (int j = 0; j < 3; j++)
	{
		double axis[3] = { 0, 0, 0 };
		axis[j] = scale / 2.0;
		vtkMath::Add(point, axis, coord);
		interpolator->Interpolate(coord, &grad1);
		vtkMath::Subtract(point, axis, coord);
		interpolator->Interpolate(coord, &grad2);
		gradient[j] = grad1 - grad2;
	}
}

void RayFeatures(vtkImageInterpolator* interpolator, const double point[3], const double normal[3], double thickness, cv::Mat& features)
{
	double coord[3];
	double intensity, gradient[3], gradmag;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 3; j++) coord[j] = point[j] + thickness*normal[j] * (i - 2);
		interpolator->Interpolate(coord, &intensity);
		ImageGradient(interpolator, coord, gradient);
		gradmag = vtkMath::Norm(gradient);
		features.at<float>(5 * i) = (float)intensity;
		features.at<float>(5 * i + 1) = (float)intensity*intensity;
		features.at<float>(5 * i + 2) = (float)gradmag;
		features.at<float>(5 * i + 3) = (float)gradmag*gradmag;
		features.at<float>(5 * i + 4) = (float)vtkMath::Dot(gradient, normal);
	}
}

int CentralizedThisContour(double center[3], double axis1[3], double axis2[3], int RingSize, double center_new[3], double* Radius, double* Thickness)
{
	vtkSmartPointer<vtkPoints> cirpoints = vtkSmartPointer<vtkPoints>::New();
	//vtkPolygon* cirPolygon = vtkPolygon::New();
	vtkSmartPointer<vtkPolyLine> cirPolyLine = vtkSmartPointer<vtkPolyLine>::New();

	double cirstep = 2.0 * M_PI / RingSize;
	for (int l = 0; l < 3; l++)
		center_new[l] = 0.0;

	for (int k = 0; k < RingSize; k++)
	{
		double ray[3], coord[3];
		for (int l = 0; l < 3; l++)
		{
			ray[l] = cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l];
			coord[l] = center[l] + Radius[k] * ray[l];
			center_new[l] = center_new[l] + coord[l];
		}
		cirpoints->InsertNextPoint(coord);
	}
	for (int l = 0; l < 3; l++)
		center_new[l] = center_new[l] / RingSize;

	cirPolyLine->GetPoints()->DeepCopy(cirpoints);
	cirPolyLine->GetPointIds()->SetNumberOfIds(RingSize);
	for (int k = 0; k < RingSize; k++)
	{
		cirPolyLine->GetPointIds()->SetId(k, k);
	}

	double tolerance = 0.01;
	double t; // Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
	double x[3]; // The coordinate of the intersection
	double pcoords[3];
	int subId;

	for (int k = 0; k < RingSize; k++)
	{
		double coord[3];
		for (int l = 0; l < 3; l++)
			coord[l] = center_new[l] + 20.0 * (cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l]);

		vtkIdType iD = cirPolyLine->IntersectWithLine(center_new, coord, tolerance, t, x, pcoords, subId);
		Radius[k] = 20.0 * t;
	}
	return 0;
}

bool PinPolyX(std::vector<double> *poly, int x, int y)
{
	bool inside = false;

	if (poly[0][y * 2] >= x && poly[0][y * 2 + 1] < x)
		inside = !inside;

	return inside;
}

//point in polygon 10.19.2014
//nVert: number of vertices in the polygon
//dVertX, dVertY: coordinates of the polygon's vertices
//x, y: coordinate of the cross section point
bool PinPoly(int nVert, std::vector<double> dVertX, std::vector<double> dVertY, double x, double y)
{
	int i, j = nVert - 1;
	bool inside = false;
	for (i = 0; i < nVert; i++)
	{
		if (dVertY[i] < y && dVertY[j] >= y || dVertY[j] < y && dVertY[i] >= y)
		{
			if (dVertX[i] + (y - dVertY[i]) / (dVertY[j] - dVertY[i]) * (dVertX[j] - dVertX[i]) < x)
			{
				inside = !inside;
			}
		}
		j = i;
	}
	return inside;
}

//mask of plaque 10.24.2014
bool PinPoly(std::vector<double> *poly, int x, int y)
{
	int i, j = poly[0].size() - 1;
	bool inside = false;
	for (i = 0; i < poly[0].size(); i++)
	{
		if (poly[1][i] < y && poly[1][j] >= y || poly[1][j] < y && poly[1][i] >= y)
		{
			if (poly[0][i] + (y - poly[1][i]) / (poly[1][j] - poly[1][i]) * (poly[0][j] - poly[0][i]) < x)
			{
				inside = !inside;
			}
		}
		j = i;
	}
	return inside;
}