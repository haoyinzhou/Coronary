#include "ImageObliqueReformat.h"

#include <vector>
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkShortArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkImageInterpolator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "common.h"
#define SAMPLING_SIZE 1.0
vtkStandardNewMacro(ImageObliqueReformat);


ImageObliqueReformat::ImageObliqueReformat( )
{
	this->SegmentId  = -1;
	this->PointId = 0;

	this->RadialSpacing  = 0.5/SAMPLING_SIZE;
	this->RadialExtent   = 20*SAMPLING_SIZE;

	this->UpdateImage = 1;

	this->SetNumberOfInputPorts( 2 );
	this->SetNumberOfOutputPorts( 7 );

	// by default process active point scalars
	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		vtkDataSetAttributes::SCALARS);
}

ImageObliqueReformat::~ImageObliqueReformat( )
{
}

void ImageObliqueReformat::SetUpdateImage(int update)
{
	if(this->UpdateImage != update) this->UpdateImage = update;
	//Don't set the modiflied flag on purpose
}

vtkImageData* ImageObliqueReformat::GetOutput()
{
	return vtkImageData::SafeDownCast(this->GetOutput(0));
}

vtkDataObject* ImageObliqueReformat::GetOutput(int port)
{
	return this->GetOutputDataObject(port);
}

int ImageObliqueReformat::ProcessRequest(vtkInformation* request,
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
int ImageObliqueReformat::FillInputPortInformation(int port, vtkInformation *info)
{
	if( port == 0 )
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
	else if( port == 1 )
		info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");

	return 1;
}

//----------------------------------------------------------------------------
int ImageObliqueReformat::FillOutputPortInformation(
	int port, vtkInformation* info)
{
	if( port == 0 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
	else if ( port < 6 )
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	else if (port == 6 ) //Mask Image
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
	return 1;
}

int ImageObliqueReformat::RequestInformation (
	vtkInformation * request,
	vtkInformationVector** inputVector,
	vtkInformationVector *outputVector)
{
	vtkInformation *inCenterlineInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outImageInfo	 = outputVector->GetInformationObject(0);

	vtkPolyData	   *inputCenterline  = vtkPolyData::SafeDownCast(inCenterlineInfo->Get(vtkDataObject::DATA_OBJECT()));

	if( SegmentId >= 0 && SegmentId < inputCenterline->GetNumberOfCells() )
	{
		int RadialSize = 2*this->RadialExtent+1;
		int outWholeExt[6] = {0, RadialSize-1, 0, RadialSize-1, 0, 0};
		outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	}
	else
	{
		int outWholeExt[6] = {0, 0, 0, 0, 0, 0};
		outImageInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),outWholeExt,6);
	}
	double outOrigin[3] = {0.0, 0.0, 0.0};
	outImageInfo->Set(vtkDataObject::ORIGIN(), outOrigin, 3);

	return Superclass::RequestInformation(request, inputVector, outputVector);
}

int ImageObliqueReformat::RequestData(
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
	vtkInformation *outImage2Info = outputVector->GetInformationObject(6); //Mask Image

	// get the input and output
	vtkImageData *inputImage	  = vtkImageData::SafeDownCast(inImageInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData  *inputCenterline = vtkPolyData::SafeDownCast(inCenterlineInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkImageData *outputImage      = vtkImageData::SafeDownCast(outImageInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *outImagePoly = vtkPolyData::SafeDownCast(outImagePolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out0Poly = vtkPolyData::SafeDownCast(out0PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out1Poly = vtkPolyData::SafeDownCast(out1PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out2Poly = vtkPolyData::SafeDownCast(out2PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *out3Poly = vtkPolyData::SafeDownCast(out3PolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkImageData *outputImage2	= vtkImageData::SafeDownCast(outImage2Info->Get(vtkDataObject::DATA_OBJECT()));

	vtkDoubleArray *clDir = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Dir"));
	vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Axis2"));
	vtkDoubleArray *clRadius = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("Radius"));
	vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *clWallThickness = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("WallThickness"));
	vtkDoubleArray *clLongiParam = vtkDoubleArray::SafeDownCast(inputCenterline->GetPointData()->GetArray("LongiParam"));
	vtkDoubleArray *clCircumParam = vtkDoubleArray::SafeDownCast(inputCenterline->GetCellData()->GetArray("CircumParam"));

	// Compute the output centerline and centerline frames
	if( SegmentId < 0 || SegmentId >= inputCenterline->GetNumberOfCells() || !clDir || !clAxis1 || !clAxis2 || !clLumenRadius || !clWallThickness || !clLongiParam || !clCircumParam )
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

	vtkIdType npts=0, *pts=NULL;
	inputCenterline->BuildCells();
	inputCenterline->GetCellPoints(SegmentId, npts, pts);
	if( PointId < 0 ) PointId=0;
	else if( PointId >= npts ) PointId = npts-1;
	vtkIdType pointId = pts[this->PointId];

	if(this->UpdateImage)
	{
		vtkImageInterpolator *interpolator = vtkImageInterpolator::New();
		interpolator->SetInterpolationModeToLinear();
		interpolator->SetOutValue(-3024.0);
		interpolator->Initialize(inputImage);

		int RadialSize = 2*this->RadialExtent+1;
		outputImage->SetExtent(0, RadialSize-1, 0, RadialSize-1, 0, 0);
		outputImage->SetOrigin(0.0, 0.0, 0.0);
		outputImage->SetSpacing(this->RadialSpacing, this->RadialSpacing, 1.0);
		outputImage->AllocateScalars(VTK_SHORT, 1);

		outputImage2->SetExtent(0, RadialSize-1, 0, RadialSize-1, 0, 0);
		outputImage2->SetOrigin(0.0, 0.0, 0.0);
		outputImage2->SetSpacing(this->RadialSpacing, this->RadialSpacing, 1.0);
		outputImage2->AllocateScalars(VTK_SHORT, 1);

		short* pixel = static_cast<short*>(outputImage->GetScalarPointer());		
		short* pixel2 = static_cast<short*>(outputImage2->GetScalarPointer());		
		double ivalue;
		double axis1[3], axis2[3];
		double coord[3], coord1[3];
		vtkIdType index;

		vtkPoints *imPoints = vtkPoints::New();
		imPoints->SetNumberOfPoints(RadialSize*RadialSize);
		vtkSmartPointer<vtkPoints> imPoints_for3DSlicer = vtkSmartPointer<vtkPoints>::New();


		outImagePoly->Allocate((RadialSize-1)*(RadialSize-1));
		vtkShortArray *imPixels = vtkShortArray::New();
		imPixels->SetNumberOfValues(imPoints->GetNumberOfPoints());
		imPixels->SetName("pointcolor");

		inputCenterline->GetPoint(pointId, coord1);
		clAxis1->GetTuple(pointId, axis1);
		clAxis2->GetTuple(pointId, axis2);

		for(vtkIdType i=-this->RadialExtent; i<=this->RadialExtent; i++)
		{
			for(vtkIdType j=-this->RadialExtent; j<=this->RadialExtent; j++)
			{
				for(int k=0; k<3; k++) coord[k] = coord1[k] + i*this->RadialSpacing*axis2[k]+ j*this->RadialSpacing*axis1[k];
				interpolator->Interpolate(coord, &ivalue);
				index = (i+this->RadialExtent)*RadialSize+j+this->RadialExtent;
				pixel[index] = short(ivalue);
				pixel2[index] = short(ivalue);
				imPoints->SetPoint(index, coord);
				imPixels->SetValue(index, short(ivalue));				
			}
		}
		
		interpolator->Delete();

		vtkIdType ids[4];
		for(vtkIdType i=-this->RadialExtent; i<this->RadialExtent; i++)
		{
			for(vtkIdType j=-this->RadialExtent; j<this->RadialExtent; j++)
			{
				ids[0] = (i+this->RadialExtent)*RadialSize+j+this->RadialExtent;
				ids[1] = ((i+1)+this->RadialExtent)*RadialSize+j+this->RadialExtent;
				ids[2] = ((i+1)+this->RadialExtent)*RadialSize+(j+1)+this->RadialExtent;
				ids[3] = (i+this->RadialExtent)*RadialSize+(j+1)+this->RadialExtent;
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
	out0Points->SetNumberOfPoints(clLumenRadius->GetNumberOfComponents());
	vtkPoints *out1Points = vtkPoints::New();
	out1Points->SetNumberOfPoints(clLumenRadius->GetNumberOfComponents());
	vtkPoints *out2Points = vtkPoints::New();
	out2Points->SetNumberOfPoints(clLumenRadius->GetNumberOfComponents());
	vtkPoints *out3Points = vtkPoints::New();
	out3Points->SetNumberOfPoints(1);

	vtkIdList *out0IdList = vtkIdList::New();
	out0IdList->SetNumberOfIds(clLumenRadius->GetNumberOfComponents()+1);
	vtkIdList *out1IdList = vtkIdList::New();
	out1IdList->SetNumberOfIds(clLumenRadius->GetNumberOfComponents()+1);
	vtkIdList *out2IdList = vtkIdList::New();
	out2IdList->SetNumberOfIds(clLumenRadius->GetNumberOfComponents()+1);
	vtkIdList *out3IdList = vtkIdList::New();
	out3IdList->SetNumberOfIds(1);

	vtkDoubleArray *out0Param = vtkDoubleArray::New();
	out0Param->SetName("Param");
	out0Param->SetNumberOfComponents(2);
	out0Param->SetNumberOfTuples(clLumenRadius->GetNumberOfComponents());
	vtkDoubleArray *out1Param = vtkDoubleArray::New();
	out1Param->SetName("Param");
	out1Param->SetNumberOfComponents(2);
	out1Param->SetNumberOfTuples(clLumenRadius->GetNumberOfComponents());
	vtkDoubleArray *out2Param = vtkDoubleArray::New();
	out2Param->SetName("Param");
	out2Param->SetNumberOfComponents(2);
	out2Param->SetNumberOfTuples(clLumenRadius->GetNumberOfComponents());

	double *lumenRadius = new double[clLumenRadius->GetNumberOfComponents()];
	double *wallThickness = new double[clWallThickness->GetNumberOfComponents()];
	double *circumparam = new double[clCircumParam->GetNumberOfComponents()];
	double circstep  = 2.0*M_PI/clLumenRadius->GetNumberOfComponents();
	double radius = clRadius->GetValue(pointId);
	clLumenRadius->GetTuple(pointId, lumenRadius);
	clWallThickness->GetTuple(pointId, wallThickness);
	clCircumParam->GetTuple(SegmentId, circumparam);
	double longiparam = clLongiParam->GetValue(pointId);

	
	for(vtkIdType i=0; i<clLumenRadius->GetNumberOfComponents(); i++)
	{
		out0Points->SetPoint(i, this->RadialExtent*outputImageSpacings[0]+radius*cos(i*circstep), this->RadialExtent*outputImageSpacings[1]+radius*sin(i*circstep), 0.0);
		out0IdList->SetId(i, i);
		out0Param->SetTuple2(i, longiparam, circumparam[i]);
		
		//lumen
		out1Points->SetPoint(i, this->RadialExtent*outputImageSpacings[0]+lumenRadius[i]*cos(i*circstep), this->RadialExtent*outputImageSpacings[1]+lumenRadius[i]*sin(i*circstep), 0.0);
		out1IdList->SetId(i, i);
		out1Param->SetTuple2(i, longiparam, circumparam[i]);
		
		//wall
		out2Points->SetPoint(i, this->RadialExtent*outputImageSpacings[0]+(lumenRadius[i]+wallThickness[i])*cos(i*circstep), this->RadialExtent*outputImageSpacings[1]+(lumenRadius[i]+wallThickness[i])*sin(i*circstep), 0.0);
		out2IdList->SetId(i, i);
		out2Param->SetTuple2(i, longiparam, circumparam[i]);			
	}	
	
	out0IdList->SetId(clLumenRadius->GetNumberOfComponents(), out0IdList->GetId(0));
	out1IdList->SetId(clLumenRadius->GetNumberOfComponents(), out1IdList->GetId(0));
	out2IdList->SetId(clLumenRadius->GetNumberOfComponents(), out2IdList->GetId(0));
	out3Points->SetPoint(0, this->RadialExtent*outputImageSpacings[0], this->RadialExtent*outputImageSpacings[1], 0.0);
	out3IdList->SetId(0, 0);

	delete[] lumenRadius;
	delete[] wallThickness;
	delete[] circumparam;
	
	out0Poly->SetPoints(out0Points); out0Points->Delete();
	out1Poly->SetPoints(out1Points); out1Points->Delete();
	out2Poly->SetPoints(out2Points); out2Points->Delete();
	out3Poly->SetPoints(out3Points); out3Points->Delete();
	out0Poly->GetPointData()->AddArray(out0Param); out0Param->Delete();
	out1Poly->GetPointData()->AddArray(out1Param); out1Param->Delete();
	out2Poly->GetPointData()->AddArray(out2Param); out2Param->Delete();
	out0Poly->Allocate(1);
	out0Poly->InsertNextCell(VTK_POLY_LINE, out0IdList); out0IdList->Delete();
	out1Poly->Allocate(1);
	out1Poly->InsertNextCell(VTK_POLY_LINE, out1IdList); out1IdList->Delete();
	out2Poly->Allocate(1);
	out2Poly->InsertNextCell(VTK_POLY_LINE, out2IdList); out2IdList->Delete();
	out3Poly->Allocate(1);
	out3Poly->InsertNextCell(VTK_VERTEX, out3IdList); out3IdList->Delete();
	
	int dims[3];
	double coord[3];
	std::vector<double> vecLumen[2];
	std::vector<double> vecWall[2];
	outputImage2->GetDimensions(dims);
	
	for(int j=0; j<out2Poly->GetNumberOfPoints(); j++)
	{
		out2Poly->GetPoint(j, coord);
		vecWall[0].push_back(coord[0]*2.0*SAMPLING_SIZE);
		vecWall[1].push_back(coord[1]*2.0*SAMPLING_SIZE);
		out1Poly->GetPoint(j, coord);
		vecLumen[0].push_back(coord[0]*2.0*SAMPLING_SIZE);
		vecLumen[1].push_back(coord[1]*2.0*SAMPLING_SIZE);
	}

	for(int j = 0; j < dims[1]; j++)
	{
		for(int i = 0; i < dims[0]; i++)
		{
			//not lumen and wall
			if(!PinPoly(vecWall, i, j))
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 1900);			
			//lumen
			else if(PinPoly(vecWall, i, j) && PinPoly(vecLumen, i, j))			
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 1900);						
		}
	}

/*
	//fixed color
	for(int j = 0; j < dims[1]; j++)
	{
		for(int i = 0; i < dims[0]; i++)
		{
			//not lumen and wall
			if(!PinPoly(vecWall, i, j))
				outputImage2->SetScalarComponentFromDouble(i, j, 0, 0, 0);
			//lumen
			else if(PinPoly(vecWall, i, j) && PinPoly(vecLumen, i, j))			
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

