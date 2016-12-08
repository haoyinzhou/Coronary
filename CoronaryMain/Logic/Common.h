#ifndef __VTK_WRAP__

#ifndef _Common_h__
#define _Common_h__

#include "vtkImageInterpolator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkImageData.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkMetaImageWriter.h"
#include "vtkImageResample.h"
#include "vtkAppendPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkCleanPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyLine.h"
#include "vtkTriangle.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangleFilter.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkButterflySubdivisionFilter.h"
#include "vtkMeshQuality.h"
#include "vtkTupleInterpolator.h"
#include "vtkCardinalSpline.h"

#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkOrImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"

#include "qprogressbar.h"
#include "QThread.h"


#include "psimpl.h"


#include "opencv2/core/core.hpp"
#include "opencv2/ml/ml.hpp"

#include "LearningImpl.h"
#include "ExtendSplineFilter.h"

using namespace cv;

namespace SmartCoronary
{
	enum LVCorLandMarkList{
		LEFT_VENTRICLE_APEX = 0,
		BASE_LEFT_END,
		BASE_RIGHT_END,
		BASE_ANTERIOR_END,
		BASE_POSTERIOR_END,
		LEFT_CORONARY_OSTIUM,
		RIGHT_CORONARY_OSTIUM,
		//////////////////////////////////////////////////////////////////////////
		NUMBER_OF_LVCOR_LANDMARKS
	};
}



void FillIntegralImage(vtkImageData* intergalImage, vtkImageData *imageData, vtkImageInterpolator* interpolator);
void IntegralImageHist(vtkImageData* intergalImage, int corner1[3], int corner2[3], double hist[6]);
void ImageFeatures(vtkImageData* intergalImage, double coord[3], cv::Mat& featureRow);
void HistSubtract(const double hist1[6], const double hist2[6], double hist3[6]);
void InsertFeature(cv::Mat& featureRow, double hist[6], int& count);
void HistNormalize(double hist[6]);
void GenerateHessianImage(vtkImageData *imageData, vtkImageData *hessianImage, vtkImageInterpolator *interpolator, QProgressBar* progressbar);
void FillSumImage(vtkImageData* sumImage, vtkImageInterpolator* interpolator);
void SumImageHist(vtkImageData* sumImage, double* sumimage, int corner1[3], int corner2[3], double& sum);
void Hessian(vtkImageData* sumImage, double* sumimage, double coord[3], double eigvalue[3], double eigvector[3][3], int cellsize);
void GetRotationMatrix(double axis[3], double angle, double rot[3][3]);
void AxisCenterline(vtkPolyData* clModel, double planenormal[3] = NULL);

void SaveVTKImage(vtkImageData *image, const char* fileName);
void SavePolyData(vtkPolyData *poly, const char* fileName);


bool DetectLandmarks_core(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator, QProgressBar* progressbar);
bool DetectCenterline_core(vtkImageData *ImageData, vtkImageData *hessianImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3], QProgressBar* progressbar);


#endif

#endif //__VTK_WRAP__