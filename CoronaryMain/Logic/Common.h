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
#include "vtkIdFilter.h"

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


class CBifurcation
{
public:
	vtkIdType CenterPID;
	double CenterCoord[3];
	std::vector< vtkIdType > VesselID;
	std::vector< int > VesselDir;
	std::vector< vtkIdType > EndfacePID;

public:
	CBifurcation();
};

struct CEndFace
{
	double x, y, z;
	std::vector<double> rx;
	std::vector<double> ry;
	std::vector<double> rz;
	std::vector<double> realrx;
	std::vector<double> realry;
	std::vector<double> realrz;
};

struct CEndFacePoint
{
	vtkIdType index[2]; // index[0] is the vessel segment id, index[1] is the point id.
	double coord[3];
	double realcoord[3];
};

struct CBifurcationTriangle
{
	std::vector<CEndFacePoint> EndFacePoint;
};


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
void RayFeatures(vtkImageInterpolator* interpolator, const double point[3], const double normal[3], double thickness, cv::Mat& features);
int CentralizedThisContour(double center[3], double axis1[3], double axis2[3], int RingSize, double center_new[3], double* Radius, double* Thickness);


void fillBifurcationTriangle(CBifurcationTriangle* t
	, vtkIdType p1_idx0, vtkIdType p1_idx1, double p1_rx, double p1_ry, double p1_rz, double p1_realrx, double p1_realry, double p1_realrz
	, vtkIdType p2_idx0, vtkIdType p2_idx1, double p2_rx, double p2_ry, double p2_rz, double p2_realrx, double p2_realry, double p2_realrz
	, vtkIdType p3_idx0, vtkIdType p3_idx1, double p3_rx, double p3_ry, double p3_rz, double p3_realrx, double p3_realry, double p3_realrz);

int findConvexPoint(int fid, int pid, int pid2, vector<CEndFace> endfaces, CBifurcationTriangle* trianglemesh_out);
int findConvexPoint_fill_big_triangle_hole(int fid1, int fid2, int fid3, int pid1, int pid2, int pid3, vector<CEndFace> endfaces, CBifurcationTriangle* trianglemesh_out);
int MergeAlgorithm(vector<CEndFace> endfaces, double bifurcationcenter[3], vector<CBifurcationTriangle>& triangles);

int smoothvtkpolydata(vtkPolyData* Poly, int iternum, int TYPE = 1);
void SaveVTKImage(vtkImageData *image, const char* fileName);
void SavePolyData(vtkPolyData *poly, const char* fileName);


bool DetectLandmarks_core(vtkImageData *imageData, Learning& learn, double landmarks[][3], vtkImageInterpolator *interpolator, QProgressBar* progressbar);
bool DetectCenterline_core(vtkImageData *ImageData, vtkImageData *hessianImage, vtkPolyData *centerlineModel, double leftOstium[3], double rightOstium[3], QProgressBar* progressbar);
bool DetectCenterlineLumenWall_core(vtkPolyData* clModel, vtkIdType selectId, vtkImageInterpolator* interpolator, Learning &learn);

bool PinPolyX(std::vector<double> *poly, int x, int y);


#endif

#endif //__VTK_WRAP__