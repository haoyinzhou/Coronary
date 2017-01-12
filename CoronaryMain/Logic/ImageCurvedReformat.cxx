#include "ImageCurvedReformat.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkShortArray.h"
#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImageInterpolator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "common.h"

#define SAMPLING_SIZE 1.0

vtkStandardNewMacro(ImageCurvedReformat);

ImageCurvedReformat::ImageCurvedReformat( )
{
	this->SegmentId = -1;

	this->TwistIndex = 0;

	this->RadialSpacing  = 0.5/SAMPLING_SIZE;
	this->RadialExtent   = 10*SAMPLING_SIZE;

	this->UpdateImage = 1;

	this->SetNumberOfInputPorts( 2 );
	this->SetNumberOfOutputPorts( 9 );

	// by default process active point scalars
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		vtkDataSetAttributes::SCALARS);
}

ImageCurvedReformat::~ImageCurvedReformat( )
{
}

void ImageCurvedReformat::SetUpdateImage(int update)
{
	if(this->UpdateImage != update) this->UpdateImage = update;
	//Don't set the modiflied flag on purpose
}

vtkImageData* ImageCurvedReformat::GetOutput()
{
	return vtkImageData::SafeDownCast(this->GetOutput(0));
}

vtkDataObject* ImageCurvedReformat::GetOutput(int port)
{
	return this->GetOutputDataObject(port);
}

int ImageCurvedReformat::ProcessRequest(vtkInformation* request,
										vtkInformationVector** inputVector,
										vtkInformationVector* outputVector)
{
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED()))
	{
		if(!this->UpdateImage)
		{
			vtkInformation* outImageInfo = outputVector->GetInformationObject(0);
			outImageInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
			vtkInformation* outImagePolyInfo = outputVector->GetInformationObject(5);
			outImagePolyInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
		}
	}

	// generate the data
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
	{
		return this->RequestData(request, inputVector, outputVector);
	}

	// execute information
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
	{
		return this->RequestInformation(request, inputVector, outputVector);
	}

	// propagate update extent
	if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
	{
		return this->RequestUpdateExtent(request, inputVector, outputVector);
	}

	return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//---------------------------------------------------------------------------
int ImageCurvedReformat::FillInputPortInformation(int port, vtkInformation *info)
{
	if( port == 0 )
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
	else if( port == 1 )
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");

	return 1;
}

//----------------------------------------------------------------------------
int ImageCurvedReformat::FillOutputPortInformation(
	int port, vtkInformation* info)
{
	if( port == 0 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
	else if ( port < 6 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	else if ( port == 6 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
	else if ( port > 6 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	return 1;
}

int ImageCurvedReformat::RequestInformation (
	vtkInformation * request,
	vtkInformationVector** inputVector,
	vtkInformationVector *outputVector)
{
	vtkInformation *inCenterlineInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outImageInfo	 = outputVector->GetInformationObject(0);

	vtkPolyData	   *inputCenterline  = vtkPolyData::SafeDownCast(inCenterlineInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkIdType npts=0, *pts=NULL;

	int RadialSize = 2*this->RadialExtent+1;
	inputCenterline->BuildCells();

	if(this->SegmentId >= 0 && this->SegmentId < inputCenterline->GetNumberOfCells())
	{
		inputCenterline->GetCellPoints(this->SegmentId, npts, pts);

		int outWholeExt[6] = {0, RadialSize-1, 0, npts-1, 0, 0};
		outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	}
	else
	{
		int outWholeExt[6] = {0, 0, 0, 0, 0, 0};
		outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	}
	//if( SegmentId >= 0 && SegmentId < inputCenterline->GetNumberOfCells() )
	//{
	//	vtkIdType npts=0, *pts=NULL;
	//	int RadialSize = 2*this->RadialExtent+1;
	//	inputCenterline->BuildCells();
	//	inputCenterline->GetCellPoints(SegmentId, npts, pts);
	//	int outWholeExt[6] = {0, RadialSize-1, 0, npts-1, 0, 0};
	//	outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	//}
	//else
	//{
	//	int outWholeExt[6] = {0, 0, 0, 0, 0, 0};
	//	outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	//}
	
	double outOrigin[3] = {0.0, 0.0, 0.0};
	outImageInfo->Set(vtkDataObject::ORIGIN(), outOrigin, 3);

	return Superclass::RequestInformation(request, inputVector, outputVector);
}

int ImageCurvedReformat::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inImageInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *inCenterlineInfo = inputVector[1]->GetInformationObject(0);

	vtkInformation *outImageInfo = outputVector->GetInformationObject(0);
	vtkInformation *outImagePolyInfo = outputVector->GetInformationObject(5);
	vtkInformation *out0PolyInfo = outputVector->GetInformationObject(1);
	vtkInformation *out1PolyInfo = outputVector->GetInformationObject(2);
	vtkInformation *out2PolyInfo = outputVector->GetInformationObject(3);
	vtkInformation *out3PolyInfo = outputVector->GetInformationObject(4);
	vtkInformation *outImage2Info = outputVector->GetInformationObject(6);
	vtkInformation *out4PolyInfo = outputVector->GetInformationObject(7);
	vtkInformation *out5PolyInfo = outputVector->GetInformationObject(8);

	// get the input and output
	vtkImageData *inputImage	  = vtkImageData::SafeDownCast(inImageInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData  *inputCenterline = vtkPolyData::SafeDownCast(inCenterlineInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkImageData *outputImage      = vtkImageData::SafeDownCast(outImageInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *outImagePoly = vtkPolyData::SafeDownCast(outImagePolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out0Poly = vtkPolyData::SafeDownCast(out0PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out1Poly = vtkPolyData::SafeDownCast(out1PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out2Poly = vtkPolyData::SafeDownCast(out2PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out3Poly = vtkPolyData::SafeDownCast(out3PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out4Poly = vtkPolyData::SafeDownCast(out4PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out5Poly = vtkPolyData::SafeDownCast(out5PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkImageData *outputImage2      = vtkImageData::SafeDownCast(outImage2Info->Get(vtkDataObject::DATA_OBJECT()));

	vtkDoubleArray *clDir = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Dir"));
	vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Axis2"));
	vtkDoubleArray *clRadius = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Radius"));
	vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *clWallThickness = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("WallThickness"));
	vtkDoubleArray *clLongiParam = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("LongiParam"));
	vtkDoubleArray *clCircumParam = vtkDoubleArray::SafeDownCast(inputCenterline->GetCellData()->GetArray("CircumParam"));


	// Compute the output centerline and centerline frames
	if( this->SegmentId < 0 || this->SegmentId >= inputCenterline->GetNumberOfCells() || !clDir || !clAxis1 || !clAxis2 || !clLumenRadius || !clWallThickness || !clLongiParam || !clCircumParam )
	{
		outputImage->SetExtent(0, 0, 0, 0, 0, 0);
		outputImage->SetOrigin(0.0, 0.0, 0.0);
		outputImage->SetSpacing(1.0, 1.0, 1.0);
		outputImage->AllocateScalars(VTK_SHORT, 1);

		outputImage2->SetExtent(0, 0, 0, 0, 0, 0);
		outputImage2->SetOrigin(0.0, 0.0, 0.0);
		outputImage2->SetSpacing(1.0, 1.0, 1.0);
		outputImage2->AllocateScalars(VTK_SHORT, 1);

		return 1;
	}
	

	vtkIdType npts = 0, *pts = NULL;
	inputCenterline->BuildCells();
	inputCenterline->GetCellPoints(this->SegmentId, npts, pts);


	double length = 0.0;
	double coord1[3], coord2[3];
		
	for(vtkIdType i=0; i<npts; i++)
	{
		inputCenterline->GetPoints()->GetPoint(pts[i], coord1);
		if(i>0)
		{
			inputCenterline->GetPoint(pts[i-1], coord2);
			length += sqrt(vtkMath::Distance2BetweenPoints(coord1, coord2));
		}
	}	

	int twistindex = this->TwistIndex%clLumenRadius->GetNumberOfComponents();
	if(twistindex < 0)
		twistindex += clLumenRadius->GetNumberOfComponents();
	int oppositeindex = (twistindex+clLumenRadius->GetNumberOfComponents()/2)%clLumenRadius->GetNumberOfComponents();
	//int extent[6];
	//outputImage->GetExtent(extent);
	//std::cout << "updateimage: " << this->UpdateImage << " | " << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] <<  " | ";

	vtkPoints *imPoints = vtkPoints::New();
	vtkSmartPointer<vtkPoints> imPoints_for3DSlicer = vtkSmartPointer<vtkPoints>::New();

	if(this->UpdateImage)
	{
		vtkImageInterpolator *interpolator = vtkImageInterpolator::New();
		interpolator->SetInterpolationModeToLinear();
		interpolator->SetOutValue(-3024.0);
		interpolator->Initialize(inputImage);

		int RadialSize = 2*this->RadialExtent+1;
		outputImage->SetExtent(0, RadialSize-1, 0, npts-1, 0, 0);
		outputImage->SetOrigin(0.0, 0.0, 0.0);
		outputImage->SetSpacing(this->RadialSpacing, length/(npts-1), 1.0);
		outputImage->AllocateScalars(VTK_SHORT, 1);

		outputImage2->SetExtent(0, RadialSize-1, 0, npts-1, 0, 0);
		outputImage2->SetOrigin(0.0, 0.0, 0.0);
		outputImage2->SetSpacing(this->RadialSpacing, length/(npts-1), 1.0);
		outputImage2->AllocateScalars(VTK_SHORT, 1);

		short* pixel = static_cast<short*>(outputImage->GetScalarPointer());
		short* pixel2 = static_cast<short*>(outputImage2->GetScalarPointer());
		double ivalue;
		double coord[3], axis1[3], axis2[3];
		vtkIdType index;

//		vtkPoints *imPoints = vtkPoints::New();
		imPoints->SetNumberOfPoints(npts * RadialSize);
		vtkShortArray *imPixels = vtkShortArray::New();
		imPixels->SetNumberOfValues(imPoints->GetNumberOfPoints());
		imPixels->SetName("pointcolor");

		double radialindex = twistindex * 2.0 * M_PI / clLumenRadius->GetNumberOfComponents();

		for(vtkIdType i = 0; i < npts; i ++)
		{
			inputCenterline->GetPoint(pts[i], coord1);
			clAxis1->GetTuple(pts[i], axis1);
			clAxis2->GetTuple(pts[i], axis2);
			for(vtkIdType j = -this->RadialExtent; j <= this->RadialExtent; j ++)
			{
				for(int k = 0; k < 3; k++)
					coord[k] = coord1[k] + j * this->RadialSpacing * (cos(radialindex) * axis1[k] + sin(radialindex) * axis2[k]);
				
				interpolator->Interpolate(coord, &ivalue);
				index = i * RadialSize + j + this->RadialExtent;
				pixel[index] = short(ivalue);
				pixel2[index] = short(ivalue);
				imPoints->SetPoint(index, coord);
				imPixels->SetValue(index, short(ivalue));
			}	
		}
		interpolator->Delete();

		outImagePoly->Allocate((npts - 1) * (RadialSize - 1));
		vtkIdType ids[4];
		for(vtkIdType i=0; i<npts-1; i++)
		{
			for(vtkIdType j=-this->RadialExtent; j<this->RadialExtent; j++)
			{
				ids[0] = i*RadialSize+j+this->RadialExtent;
				ids[1] = (i+1)*RadialSize+j+this->RadialExtent;
				ids[2] = (i+1)*RadialSize+(j+1)+this->RadialExtent;
				ids[3] = i*RadialSize+(j+1)+this->RadialExtent;
				outImagePoly->InsertNextCell(VTK_QUAD, 4, ids);
			}
		}

		for (int i = 0; i < imPoints->GetNumberOfPoints(); i++)
		{
			double coordi[3];
			imPoints->GetPoint(i, coordi);
			coordi[0] = -coordi[0];
			coordi[1] = -coordi[1];
			imPoints_for3DSlicer->InsertNextPoint(coordi);
		}

		outImagePoly->SetPoints(imPoints_for3DSlicer); imPoints->Delete();
		outImagePoly->GetPointData()->AddArray(imPixels); imPixels->Delete();
	}
	double outputImageSpacings[3];
	outputImage->GetSpacing(outputImageSpacings);

	vtkPoints *out0Points = vtkPoints::New();
	out0Points->SetNumberOfPoints(2*npts);
	vtkPoints *out1Points = vtkPoints::New();
	out1Points->SetNumberOfPoints(2*npts);
	vtkPoints *out2Points = vtkPoints::New();
	out2Points->SetNumberOfPoints(2*npts);
	vtkPoints *out3Points = vtkPoints::New();
	out3Points->SetNumberOfPoints(npts);
	vtkPoints *out4Points = vtkPoints::New();
	out4Points->SetNumberOfPoints(2*npts);
	vtkPoints *out5Points = vtkPoints::New();
	out5Points->SetNumberOfPoints(2*npts);

	vtkIdList *out0IdList1 = vtkIdList::New();
	out0IdList1->SetNumberOfIds(npts);
	vtkIdList *out0IdList2 = vtkIdList::New();
	out0IdList2->SetNumberOfIds(npts);
	vtkIdList *out1IdList1 = vtkIdList::New();
	out1IdList1->SetNumberOfIds(npts);
	vtkIdList *out1IdList2 = vtkIdList::New();
	out1IdList2->SetNumberOfIds(npts);
	vtkIdList *out2IdList1 = vtkIdList::New();
	out2IdList1->SetNumberOfIds(npts);
	vtkIdList *out2IdList2 = vtkIdList::New();
	out2IdList2->SetNumberOfIds(npts);
	vtkIdList *out3IdList = vtkIdList::New();
	out3IdList->SetNumberOfIds(npts);
	vtkIdList *out4IdList1 = vtkIdList::New();
	out4IdList1->SetNumberOfIds(npts);
	vtkIdList *out4IdList2 = vtkIdList::New();
	out4IdList2->SetNumberOfIds(npts);
	vtkIdList *out5IdList1 = vtkIdList::New();
	out5IdList1->SetNumberOfIds(npts);
	vtkIdList *out5IdList2 = vtkIdList::New();
	out5IdList2->SetNumberOfIds(npts);

	vtkDoubleArray *out0Param = vtkDoubleArray::New();
	out0Param->SetName("Param");
	out0Param->SetNumberOfComponents(2);
	out0Param->SetNumberOfTuples(2*npts);
	vtkDoubleArray *out1Param = vtkDoubleArray::New();
	out1Param->SetName("Param");
	out1Param->SetNumberOfComponents(2);
	out1Param->SetNumberOfTuples(2*npts);
	vtkDoubleArray *out2Param = vtkDoubleArray::New();
	out2Param->SetName("Param");
	out2Param->SetNumberOfComponents(2);
	out2Param->SetNumberOfTuples(2*npts);
	vtkDoubleArray *out4Param = vtkDoubleArray::New();
	out4Param->SetName("Param");
	out4Param->SetNumberOfComponents(2);
	out4Param->SetNumberOfTuples(2*npts);
	vtkDoubleArray *out5Param = vtkDoubleArray::New();
	out5Param->SetName("Param");
	out5Param->SetNumberOfComponents(2);
	out5Param->SetNumberOfTuples(2*npts);

	double *lumenRadius = new double[clLumenRadius->GetNumberOfComponents()];
	double *wallThickness = new double[clWallThickness->GetNumberOfComponents()];
	double *circumparam = new double[clCircumParam->GetNumberOfComponents()];	
	clCircumParam->GetTuple(this->SegmentId, circumparam);
	
	for(vtkIdType i = 0; i < npts; i ++)
	{
		double radius = clRadius->GetValue(pts[i]);
		clLumenRadius->GetTuple(pts[i], lumenRadius);
		clWallThickness->GetTuple(pts[i], wallThickness);
		double longiparam = clLongiParam->GetValue(pts[i]);
		
		out0Points->SetPoint(2*i, this->RadialExtent*outputImageSpacings[0]+radius, i*outputImageSpacings[1], 0.0);
		out0IdList1->SetId(i,2*i);
		out0Param->SetTuple2(2*i, longiparam, circumparam[twistindex]);

		out0Points->SetPoint(2*i+1, this->RadialExtent*outputImageSpacings[0]-radius, i*outputImageSpacings[1], 0.0);
		out0IdList2->SetId(i,2*i+1);
		out0Param->SetTuple2(2*i+1, longiparam, circumparam[oppositeindex]);

		out1Points->SetPoint(2*i, this->RadialExtent*outputImageSpacings[0]+lumenRadius[twistindex], i*outputImageSpacings[1], 0.0);
		out1IdList1->SetId(i, 2*i);
		out1Param->SetTuple2(2*i, longiparam, circumparam[twistindex]);
		out1Points->SetPoint(2*i+1, this->RadialExtent*outputImageSpacings[0]-lumenRadius[oppositeindex], i*outputImageSpacings[1], 0.0);
		out1IdList2->SetId(i, 2*i+1);
		out1Param->SetTuple2(2*i+1, longiparam, circumparam[oppositeindex]);

		out2Points->SetPoint(2*i, this->RadialExtent*outputImageSpacings[0]+(lumenRadius[twistindex]+wallThickness[twistindex]), i*outputImageSpacings[1], 0.0);
		out2IdList1->SetId(i, 2*i);
		out2Param->SetTuple2(2*i, longiparam, circumparam[twistindex]);
		out2Points->SetPoint(2*i+1, this->RadialExtent*outputImageSpacings[0]-(lumenRadius[oppositeindex]+wallThickness[oppositeindex]), i*outputImageSpacings[1], 0.0);
		out2IdList2->SetId(i, 2*i+1);
		out2Param->SetTuple2(2*i+1, longiparam, circumparam[oppositeindex]);

		out3Points->SetPoint(i, this->RadialExtent*outputImageSpacings[0], i*outputImageSpacings[1], 0.0);
		out3IdList->SetId(i, i);
	}

//	int RadialSize = 2*this->RadialExtent+1;
	double radialindex = twistindex*2.0*M_PI/clLumenRadius->GetNumberOfComponents();	
	double coord[3], axis1[3], axis2[3];

	for(vtkIdType i = 0; i < npts; i ++)
	{
//		double radius = clRadius->GetValue(pts[i]);
		clLumenRadius->GetTuple(pts[i], lumenRadius);
		clWallThickness->GetTuple(pts[i], wallThickness);
		double longiparam = clLongiParam->GetValue(pts[i]);
		double cstep = 1.0/clLumenRadius->GetNumberOfComponents();
		double cirstep = 2.0*M_PI*cstep;

		inputCenterline->GetPoint(pts[i], coord1);
		clAxis1->GetTuple(pts[i], axis1);
		clAxis2->GetTuple(pts[i], axis2);

		for(int k=0; k<3; k++) coord[k] = coord1[k] + lumenRadius[twistindex]*(cos(twistindex*cirstep)*axis1[k]+sin(twistindex*cirstep)*axis2[k]) ;
		out4Points->SetPoint(2*i, coord[0], coord[1], coord[2]);
		out4IdList1->SetId(i,2*i);
		out4Param->SetTuple2(2*i, longiparam, circumparam[twistindex]);
		for(int k=0; k<3; k++) coord[k] = coord1[k] + lumenRadius[oppositeindex]*(cos(oppositeindex*cirstep)*axis1[k]+sin(oppositeindex*cirstep)*axis2[k]) ;
		out4Points->SetPoint(2*i+1, coord[0], coord[1], coord[2]);
		out4IdList2->SetId(i,2*i+1);
		out4Param->SetTuple2(2*i+1, longiparam, circumparam[oppositeindex]);

		for(int k=0; k<3; k++) coord[k] = coord1[k] + (wallThickness[twistindex]+lumenRadius[twistindex])*(cos(twistindex*cirstep)*axis1[k]+sin(twistindex*cirstep)*axis2[k]) ;
		out5Points->SetPoint(2*i, coord[0], coord[1], coord[2]);
		out5IdList1->SetId(i,2*i);
		out5Param->SetTuple2(2*i, longiparam, circumparam[twistindex]);
		for(int k=0; k<3; k++) coord[k] = coord1[k] + (wallThickness[oppositeindex]+lumenRadius[oppositeindex])*(cos(oppositeindex*cirstep)*axis1[k]+sin(oppositeindex*cirstep)*axis2[k]) ;
		out5Points->SetPoint(2*i+1, coord[0], coord[1], coord[2]);
		out5IdList2->SetId(i,2*i+1);
		out5Param->SetTuple2(2*i+1, longiparam, circumparam[oppositeindex]);
		
	}

	delete[] lumenRadius;
	delete[] wallThickness;
	delete[] circumparam;
	
	out0Poly->SetPoints(out0Points); out0Points->Delete();
	out1Poly->SetPoints(out1Points); out1Points->Delete();
	out2Poly->SetPoints(out2Points); out2Points->Delete();
	out3Poly->SetPoints(out3Points); out3Points->Delete();
	out4Poly->SetPoints(out4Points); out4Points->Delete();
	out5Poly->SetPoints(out5Points); out5Points->Delete();
	out0Poly->GetPointData()->AddArray(out0Param); out0Param->Delete();
	out1Poly->GetPointData()->AddArray(out1Param); out1Param->Delete();
	out2Poly->GetPointData()->AddArray(out2Param); out2Param->Delete();
	out4Poly->GetPointData()->AddArray(out4Param); out4Param->Delete();
	out5Poly->GetPointData()->AddArray(out5Param); out5Param->Delete();
	out0Poly->Allocate(2);
	out0Poly->InsertNextCell(VTK_POLY_LINE, out0IdList1); out0IdList1->Delete();
	out0Poly->InsertNextCell(VTK_POLY_LINE, out0IdList2); out0IdList2->Delete();
	out1Poly->Allocate(2);
	out1Poly->InsertNextCell(VTK_POLY_LINE, out1IdList1); out1IdList1->Delete();
	out1Poly->InsertNextCell(VTK_POLY_LINE, out1IdList2); out1IdList2->Delete();
	out2Poly->Allocate(2);
	out2Poly->InsertNextCell(VTK_POLY_LINE, out2IdList1); out2IdList1->Delete();
	out2Poly->InsertNextCell(VTK_POLY_LINE, out2IdList2); out2IdList2->Delete();
	out3Poly->Allocate(1);
	out3Poly->InsertNextCell(VTK_POLY_LINE, out3IdList); out3IdList->Delete();
	out4Poly->Allocate(2);
	out4Poly->InsertNextCell(VTK_POLY_LINE, out4IdList1); out4IdList1->Delete();
	out4Poly->InsertNextCell(VTK_POLY_LINE, out4IdList2); out4IdList2->Delete();
	out5Poly->Allocate(2);
	out5Poly->InsertNextCell(VTK_POLY_LINE, out5IdList1); out5IdList1->Delete();
	out5Poly->InsertNextCell(VTK_POLY_LINE, out5IdList2); out5IdList2->Delete();
	

	int dims[3];
	double coord3[3];
	std::vector<double> vecLumen[2];
	std::vector<double> vecWall[2];
	outputImage2->GetDimensions(dims);

	for(int j=0; j<out2Poly->GetNumberOfPoints(); j++)
	{
		out2Poly->GetPoint(j, coord3);
		vecWall[0].push_back(coord3[0]*2.0*SAMPLING_SIZE);
		vecWall[1].push_back(coord3[1]*2.0*SAMPLING_SIZE);
		out1Poly->GetPoint(j, coord3);
		vecLumen[0].push_back(coord3[0]*2.0*SAMPLING_SIZE);
		vecLumen[1].push_back(coord3[1]*2.0*SAMPLING_SIZE);
	}

	for(int j = 0; j < dims[1]; j++)
	{
		for(int i = 0; i < dims[0]; i++)
		{
			//not lumen and wall
			if(!PinPolyX(vecWall, i, j))
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 1900);
			//lumen
			else if(PinPolyX(vecWall, i, j) && PinPolyX(vecLumen, i, j))			
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 1900);			
		}
	}

/*
	for(int j = 0; j < dims[1]; j++)
	{
		for(int i = 0; i < dims[0]; i++)
		{
			//not lumen and wall
			if(!PinPolyX(vecWall, i, j))
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 0);
			//lumen
			else if(PinPolyX(vecWall, i, j) && PinPolyX(vecLumen, i, j))			
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 5);			
			//necrotic
			else if( -30 < outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) &&
				outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) < 30)
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 1);
			//fibrous fatty
			else if( 31 < outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) &&
				outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) < 130)
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 2);
			//fibrous
			else if( 131 < outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) &&
				outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) < 350)
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 3);
			//calcified
			else if( 350 < outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) &&
				outputImage2->GetScalarComponentAsDouble(i, j, 0, 0) < 2048)
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 4);
			//undefined
			else
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 6);
		}
	}
*/

	return 1;
}

