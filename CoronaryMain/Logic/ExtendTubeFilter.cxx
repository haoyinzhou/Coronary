#include "ExtendTubeFilter.h"

#include <map>
#include "common.h"


vtkStandardNewMacro(ExtendTubeFilter);

void ExtendTubeFilter::GetRotationMatrix(double axis[3], double angle, double rot[3][3])
{
	double c = cos(angle);
	double s = sin(angle);
	rot[0][0] = c + axis[0]*axis[0]*(1.0-c);
	rot[0][1] = axis[0]*axis[1]*(1.0-c) - s*axis[2];
	rot[0][2] = axis[0]*axis[2]*(1.0-c) + s*axis[1];
	rot[1][0] = axis[0]*axis[1]*(1.0-c) + s*axis[2];
	rot[1][1] = c + axis[1]*axis[1]*(1.0-c);
	rot[1][2] = axis[1]*axis[2]*(1.0-c) - s*axis[0];
	rot[2][0] = axis[0]*axis[2]*(1.0-c) - s*axis[1];
	rot[2][1] = axis[1]*axis[2]*(1.0-c) + s*axis[0];
	rot[2][2] = c + axis[2]*axis[2]*(1.0-c);
}

ExtendTubeFilter::ExtendTubeFilter()
{
	this->Capping = 0;
	this->LongitudinalRefineSteps = 3;//2;
	this->CircumferentialRefineSteps = 0;//1;
	this->LongitudinalResampleSteps = 1;
	this->CircumferentialResampleSteps = 1; // do not change it
	this->RadiusScale = 1.0;

	this->UpdateSegment = -1;
	this->firstSegment = true;

	this->SetNumberOfInputPorts( 1 );
	this->SetNumberOfOutputPorts( 5 );

	this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Radius");
	
	this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"LumenRadius");

	this->SetInputArrayToProcess(2,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"WallThickness");

	this->SetInputArrayToProcess(3,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Dir");

	this->SetInputArrayToProcess(4,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Axis1");

	this->SetInputArrayToProcess(5,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
		"Axis2");

	out0cache = vtkPolyData::New();
}

ExtendTubeFilter::~ExtendTubeFilter()
{
	out0cache->Delete();
}

void ExtendTubeFilter::SetUpdateSegment(vtkIdType update)
{
	if(this->UpdateSegment != update) this->UpdateSegment = update;
	//Don't set the modiflied flag on purpose
}

int ExtendTubeFilter::ProcessRequest(vtkInformation* request,
                                         vtkInformationVector** inputVector,
                                         vtkInformationVector* outputVector)
{
	if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED()))
	{
		if(this->UpdateSegment >= 0)
		{
			for(int i = 1; i < 5; i ++)
			{
				vtkInformation* outInfo = outputVector->GetInformationObject(i);
				outInfo->Set(vtkDemandDrivenPipeline::DATA_NOT_GENERATED(), 1);
			}
		}
	}

  // generate the data
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->RequestData(request, inputVector, outputVector);
    }

  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
    }

  // execute information
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return this->RequestInformation(request, inputVector, outputVector);
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int ExtendTubeFilter::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);	
	vtkInformation *out0Info = outputVector->GetInformationObject(0);
	vtkInformation *out1Info = outputVector->GetInformationObject(1);
	vtkInformation *out2Info = outputVector->GetInformationObject(2);
	vtkInformation *out3Info = outputVector->GetInformationObject(3);
	vtkInformation *out4Info = outputVector->GetInformationObject(4);
	

	// get the input and output
	vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));	
	vtkPolyData *output0 = vtkPolyData::SafeDownCast(out0Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output1 = vtkPolyData::SafeDownCast(out1Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output2 = vtkPolyData::SafeDownCast(out2Info->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output3 = vtkPolyData::SafeDownCast(out3Info->Get(vtkDataObject::DATA_OBJECT()));
//	vtkPolyData *output4 = vtkPolyData::SafeDownCast(out4Info->Get(vtkDataObject::DATA_OBJECT()));

	vtkDoubleArray *clRadius=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(0,inputVector));
	vtkDoubleArray *clLumenRadius=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(1,inputVector));
	vtkDoubleArray *clWallThickness=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(2,inputVector));
	vtkDoubleArray *clDir=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(3,inputVector));
	vtkDoubleArray *clAxis1=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(4,inputVector));
	vtkDoubleArray *clAxis2=vtkDoubleArray::SafeDownCast(this->GetInputArrayToProcess(5,inputVector));
	if(!clRadius)
	{
		vtkDebugMacro("Could not find clRadius.");
		return 1;
	}
	if(!clLumenRadius)
	{
		vtkDebugMacro("Could not find clLumenRadius.");
		return 1;
	}
	if(!clWallThickness)
	{
		vtkDebugMacro("Could not find clWallThickness.");
		return 1;
	}
	if( clLumenRadius->GetNumberOfComponents() != clWallThickness->GetNumberOfComponents() || clLumenRadius->GetNumberOfComponents()<3 )
	{
		vtkDebugMacro("The number of components in clLumenRadius or clWallThickness is incorrect. ");
		return 1;
	}
	if(!clDir)
	{
		vtkDebugMacro("Could not find clDir.");
		return 1;
	}
	if(!clAxis1)
	{
		vtkDebugMacro("Could not find clAxis1.");
		return 1;
	}
	if(!clAxis2)
	{
		vtkDebugMacro("Could not find clAxis2.");
		return 1;
	}

	vtkIdType inCellId;
	vtkIdType npts=0, *pts=NULL;

	vtkPoints *out0Points = vtkPoints::New();
	vtkCellArray *out0Lines = vtkCellArray::New();
	vtkDoubleArray *out0Radius = vtkDoubleArray::New();
	out0Radius->SetName("Radius");
	out0Radius->SetNumberOfComponents(1);
	vtkDoubleArray *out0LumenRadius = vtkDoubleArray::New();
	out0LumenRadius->SetName("LumenRadius");
	out0LumenRadius->SetNumberOfComponents(clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));
	vtkDoubleArray *out0WallThickness = vtkDoubleArray::New();
	out0WallThickness->SetName("WallThickness");
	out0WallThickness->SetNumberOfComponents(clWallThickness->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));
	vtkDoubleArray *out0Dir = vtkDoubleArray::New();
	out0Dir->SetName("Dir");
	out0Dir->SetNumberOfComponents(3);
	vtkDoubleArray *out0Axis1 = vtkDoubleArray::New();
	out0Axis1->SetName("Axis1");
	out0Axis1->SetNumberOfComponents(3);
	vtkDoubleArray *out0Axis2 = vtkDoubleArray::New();
	out0Axis2->SetName("Axis2");
	out0Axis2->SetNumberOfComponents(3);
	vtkDoubleArray *out0LongiParam = vtkDoubleArray::New();
	out0LongiParam->SetName("LongiParam");
	out0LongiParam->SetNumberOfComponents(1);
	vtkDoubleArray *out0CircumParam = vtkDoubleArray::New();
	out0CircumParam->SetName("CircumParam");
	out0CircumParam->SetNumberOfComponents(out0LumenRadius->GetNumberOfComponents());

	vtkCardinalSpline *spline = vtkCardinalSpline::New();

	//Cell Data for output 0
	output0->GetCellData()->CopyAllocate(input->GetCellData());
	vtkPoints*		inPoints = input->GetPoints();
	vtkCellArray* inLines = input->GetLines();

	double *radii = new double[clLumenRadius->GetNumberOfComponents()];
	double *thickness = new double[clWallThickness->GetNumberOfComponents()];
	double *refineradii = new double[clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double *refineradii_2 = new double[clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double *refinethickness = new double[clWallThickness->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double radius, coord[3], center[3];
	double sdir[3], edir[3], saxis1[3], eaxis1[3], saxis2[3], eaxis2[3];
	double dir[3], axis1[3], axis2[3];
	double longiparam;
	double *circumparam = new double[clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1)];
	double rot[3][3], rotaxis[3], rotangle;
	double cstep = 1.0/out0LumenRadius->GetNumberOfComponents();
	double cirstep  = 2.0*M_PI*cstep;
	double circumstep = clLumenRadius->GetNumberOfComponents()*cstep;

	for(inCellId=0, inLines->InitTraversal(); inLines->GetNextCell(npts,pts); inCellId++)
	{
		if( this->UpdateSegment >= 0 && inCellId != this->UpdateSegment ) continue;

		vtkIdList *idlist = vtkIdList::New();

		vtkTupleInterpolator* pointInterpolator = vtkTupleInterpolator::New();
		pointInterpolator->SetNumberOfComponents(3);
		vtkTupleInterpolator* radiusInterpolator = vtkTupleInterpolator::New();
		radiusInterpolator->SetNumberOfComponents(1);
		vtkTupleInterpolator* lumenRadiusInterpolator = vtkTupleInterpolator::New();
		lumenRadiusInterpolator->SetNumberOfComponents(clLumenRadius->GetNumberOfComponents());
		vtkTupleInterpolator* wallThicknessInterpolator = vtkTupleInterpolator::New();
		wallThicknessInterpolator->SetNumberOfComponents(clWallThickness->GetNumberOfComponents());
		for(vtkIdType j=0; j<npts; j++)
		{
			inPoints->GetPoint(pts[j], coord);
			pointInterpolator->AddTuple(j, coord);
			radius = clRadius->GetValue(pts[j]);
			radiusInterpolator->AddTuple(j, &radius);
			clLumenRadius->GetTuple(pts[j], radii);
			lumenRadiusInterpolator->AddTuple(j, radii);
			clWallThickness->GetTuple(pts[j], thickness);
			wallThicknessInterpolator->AddTuple(j, thickness);
		}

		for(vtkIdType j=0; j<npts; j++)
		{
			//out0Points, out0Radius, out0LumenRadius, out0WallThickness, out0LongiParam
			pointInterpolator->InterpolateTuple(j, coord);
			idlist->InsertNextId(out0Points->InsertNextPoint(coord));
			radiusInterpolator->InterpolateTuple(j, &radius);
			out0Radius->InsertNextValue(radius);
			lumenRadiusInterpolator->InterpolateTuple(j, radii);
			InterpolateRefine(spline, radii, clLumenRadius->GetNumberOfComponents(), refineradii, this->CircumferentialRefineSteps);
			out0LumenRadius->InsertNextTuple(refineradii);
			wallThicknessInterpolator->InterpolateTuple(j, thickness);
			InterpolateRefine(spline, thickness, clWallThickness->GetNumberOfComponents(), refinethickness, this->CircumferentialRefineSteps);
			out0WallThickness->InsertNextTuple(refinethickness);
			out0LongiParam->InsertNextValue(j);

			//out0Dir, out0Axis1, out0Axis2
			if(j==0)
			{
				if(npts==2)
				{
					inPoints->GetPoint(pts[j], coord);
					inPoints->GetPoint(pts[j+1], sdir);
					vtkMath::Subtract(sdir, coord, sdir);
					vtkMath::Normalize(sdir);
					vtkMath::Perpendiculars(sdir, saxis1, saxis2, 0.0);
				}
				else
				{
					clDir->GetTuple(pts[j+1], sdir);
					clAxis1->GetTuple(pts[j+1], saxis1);
					clAxis2->GetTuple(pts[j+1], saxis2);
				}
				std::copy(sdir, sdir+3, edir);
				std::copy(saxis1, saxis1+3, eaxis1);
				std::copy(saxis2, saxis2+3, eaxis2);
			}
			else
			{
				std::copy(edir, edir+3, sdir);
				std::copy(eaxis1, eaxis1+3, saxis1);
				std::copy(eaxis2, eaxis2+3, saxis2);
				if(j<npts-2)
				{
					clDir->GetTuple(pts[j+1], edir);
					clAxis1->GetTuple(pts[j+1], eaxis1);
					clAxis2->GetTuple(pts[j+1], eaxis2);
				}
				else
				{
					std::copy(sdir, sdir+3, edir);
					std::copy(saxis1, saxis1+3, eaxis1);
					std::copy(saxis2, saxis2+3, eaxis2);
				}
			}
			out0Dir->InsertNextTuple(sdir);
			out0Axis1->InsertNextTuple(saxis1);
			out0Axis2->InsertNextTuple(saxis2);

			if(j==npts-1) break;

			if(this->LongitudinalRefineSteps>0)
			{
				if( edir[0] == sdir[0] && edir[1] == sdir[1] && edir[2] == sdir[2] )
				{
					rotaxis[0] = 0.0; rotaxis[1] = 0.0; rotaxis[2] = 0.0;
					rotangle = 0.0;
				}
				else
				{
					vtkMath::Cross(sdir, edir, rotaxis);
					vtkMath::Normalize(rotaxis);
					rotangle = acos(vtkMath::Dot(sdir, edir));
				}

				double lrs = 1.0/(this->LongitudinalRefineSteps+1);
				for(int k=1; k<=this->LongitudinalRefineSteps; k++)
				{
					double s = j+k*lrs;
					pointInterpolator->InterpolateTuple(s, coord);
					idlist->InsertNextId(out0Points->InsertNextPoint(coord));
					radiusInterpolator->InterpolateTuple(s, &radius);
					out0Radius->InsertNextValue(radius);
					lumenRadiusInterpolator->InterpolateTuple(s, radii);
					InterpolateRefine(spline, radii, clLumenRadius->GetNumberOfComponents(), refineradii, this->CircumferentialRefineSteps);
					out0LumenRadius->InsertNextTuple(refineradii);
					wallThicknessInterpolator->InterpolateTuple(s, thickness);
					InterpolateRefine(spline, thickness, clWallThickness->GetNumberOfComponents(), refinethickness, this->CircumferentialRefineSteps);
					out0WallThickness->InsertNextTuple(refinethickness);
					out0LongiParam->InsertNextValue(s);

					double angle = k*lrs*rotangle;
					GetRotationMatrix(rotaxis, angle, rot);
					vtkMath::Multiply3x3(rot, saxis1, axis1);
					vtkMath::Multiply3x3(rot, saxis2, axis2);
					vtkMath::Multiply3x3(rot, sdir, dir);
					vtkMath::Normalize(axis1);
					vtkMath::Normalize(axis2);
					vtkMath::Normalize(dir);
					out0Dir->InsertNextTuple(dir);
					out0Axis1->InsertNextTuple(axis1);
					out0Axis2->InsertNextTuple(axis2);
				}
			}
		}
				
		pointInterpolator->Delete();
		radiusInterpolator->Delete();
		lumenRadiusInterpolator->Delete();
		wallThicknessInterpolator->Delete();

		vtkIdType outcellId = out0Lines->InsertNextCell(idlist);
		idlist->Delete();
		output0->GetCellData()->CopyData(input->GetCellData(),inCellId,outcellId);

		for(int k=0; k<out0CircumParam->GetNumberOfComponents(); k++)
		{
			circumparam[k] = k*circumstep;
		}
		out0CircumParam->InsertNextTuple(circumparam);
	}
	delete[] radii;
	delete[] thickness;
	spline->Delete();


	if( this->UpdateSegment >= 0 )
	{
		//Assume the topology would not change
		vtkPoints* cachePoints = out0cache->GetPoints();
		vtkDoubleArray *cacheRadius = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Radius"));
		vtkDoubleArray *cacheLumenRadius = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("LumenRadius"));
		vtkDoubleArray *cacheWallThickness = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("WallThickness"));
		vtkDoubleArray *cacheDir = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Dir"));
		vtkDoubleArray *cacheAxis1 = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Axis1"));
		vtkDoubleArray *cacheAxis2 = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("Axis2"));
		vtkDoubleArray *cacheLongiParam = vtkDoubleArray::SafeDownCast(out0cache->GetPointData()->GetArray("LongiParam"));
		vtkDoubleArray *cacheCircumParam = vtkDoubleArray::SafeDownCast(out0cache->GetCellData()->GetArray("CircumParam"));
		inLines = out0cache->GetLines();
		for(inCellId=0, inLines->InitTraversal(); inLines->GetNextCell(npts,pts); inCellId++)
		{
			if(inCellId == this->UpdateSegment)
			{
				for(vtkIdType j=0; j<npts; j++)
				{
					cachePoints->SetPoint(pts[j], out0Points->GetPoint(j));
					cacheRadius->SetValue(pts[j], out0Radius->GetValue(j));
					cacheLumenRadius->SetTuple(pts[j], out0LumenRadius->GetTuple(j));
					cacheWallThickness->SetTuple(pts[j], out0WallThickness->GetTuple(j));
					cacheDir->SetTuple(pts[j], out0Dir->GetTuple(j));
					cacheAxis1->SetTuple(pts[j], out0Axis1->GetTuple(j));
					cacheAxis2->SetTuple(pts[j], out0Axis2->GetTuple(j));
					cacheLongiParam->SetValue(pts[j], out0LongiParam->GetValue(j));
				}
				cacheCircumParam->SetTuple(inCellId, out0CircumParam->GetTuple(0));
				break;
			}
		}
		output0->DeepCopy(out0cache);
		
	}
	else
	{
		output0->SetPoints(out0Points);
		output0->SetLines(out0Lines);
		output0->GetPointData()->AddArray(out0Radius);
		output0->GetPointData()->AddArray(out0LumenRadius);
		output0->GetPointData()->AddArray(out0WallThickness);
		output0->GetPointData()->AddArray(out0Dir);
		output0->GetPointData()->AddArray(out0Axis1);
		output0->GetPointData()->AddArray(out0Axis2);
		output0->GetPointData()->AddArray(out0LongiParam);
		output0->GetCellData()->AddArray(out0CircumParam);
	}
	out0cache->DeepCopy(output0);

	if( this->UpdateSegment < 0 )
	{

	// by Haoyin Zhou
	printf("ExtendTubeFilter begin!\n");

	int NumofInsectParts = 0;
	CInsectPart* InsertPart = new CInsectPart[200];
	CInsectPart InsertPartBeg, InsertPartEnd;
	bool* EndFaceisfind = new bool[200];

	vtkPolyData* polydataCopy = vtkPolyData::New();
	polydataCopy->DeepCopy(output0);

	vtkCellArray* out0Lines2 = polydataCopy->GetLines();
	vtkIdType outCellId = 0, outCellId2 = 0;
	vtkIdType npts2 = 0, *pts2 = NULL;
	bool BegAlreadyIn = false, EndAlreadyIn = false;
	
	double coordbeg1[3], coordend1[3], coordbeg2[3], coordend2[3];
	for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
	{
		BegAlreadyIn = false;
		EndAlreadyIn = false;

		InsertPartBeg.NumofSeg = 1;
		InsertPartBeg.ts = pts[0];
		out0Points->GetPoint(InsertPartBeg.ts, InsertPartBeg.coord);
		InsertPartBeg.SegIDs[0] = outCellId;
		InsertPartBeg.SegDir[0] = 1;

		InsertPartEnd.NumofSeg = 1;
		InsertPartEnd.ts = pts[npts - 1];
		out0Points->GetPoint(InsertPartEnd.ts, InsertPartEnd.coord);
		InsertPartEnd.SegIDs[0] = outCellId;
		InsertPartEnd.SegDir[0] = -1;

		for (int i = 0; i < NumofInsectParts; i ++)
		{
			if (abs(InsertPartBeg.coord[0] - InsertPart[i].coord[0]) < 1e-6
			 && abs(InsertPartBeg.coord[1] - InsertPart[i].coord[1]) < 1e-6
			 && abs(InsertPartBeg.coord[2] - InsertPart[i].coord[2]) < 1e-6)
			{
				BegAlreadyIn = true;
				break;
			}
		}
		for (int i = 0; i < NumofInsectParts; i ++)
		{
			if (abs(InsertPartEnd.coord[0] - InsertPart[i].coord[0]) < 1e-6
			 && abs(InsertPartEnd.coord[1] - InsertPart[i].coord[1]) < 1e-6
			 && abs(InsertPartEnd.coord[2] - InsertPart[i].coord[2]) < 1e-6)
			{
				EndAlreadyIn = true;
				break;
			}
		}

		if (BegAlreadyIn == true && EndAlreadyIn == true)
		{
			continue;
		}

		out0Points->GetPoint(pts[0], coordbeg1);
		out0Points->GetPoint(pts[npts - 1], coordend1);

		for(outCellId2 = 0, out0Lines2->InitTraversal(); out0Lines2->GetNextCell(npts2, pts2); outCellId2 ++)
		{
			if (outCellId2 == outCellId)
				continue;

			out0Points->GetPoint(pts2[0], coordbeg2);
			out0Points->GetPoint(pts2[npts2-1], coordend2);

			if (abs(coordbeg1[0] - coordbeg2[0]) < 1e-6 && abs(coordbeg1[1] - coordbeg2[1]) < 1e-6 && abs(coordbeg1[2] - coordbeg2[2]) < 1e-6)
			{
				InsertPartBeg.SegIDs[InsertPartBeg.NumofSeg] = outCellId2;
				InsertPartBeg.SegDir[InsertPartBeg.NumofSeg] = 1;
				InsertPartBeg.NumofSeg ++;
			}
			else if (abs(coordbeg1[0] - coordend2[0]) < 1e-6 && abs(coordbeg1[1] - coordend2[1]) < 1e-6 && abs(coordbeg1[2] - coordend2[2]) < 1e-6)
			{
				InsertPartBeg.SegIDs[InsertPartBeg.NumofSeg] = outCellId2;
				InsertPartBeg.SegDir[InsertPartBeg.NumofSeg] = -1;
				InsertPartBeg.NumofSeg ++;
			}

			if (abs(coordend1[0] - coordbeg2[0]) < 1e-6 && abs(coordend1[1] - coordbeg2[1]) < 1e-6 && abs(coordend1[2] - coordbeg2[2]) < 1e-6)
			{
				InsertPartEnd.SegIDs[InsertPartEnd.NumofSeg] = outCellId2;
				InsertPartEnd.SegDir[InsertPartEnd.NumofSeg] = 1;
				InsertPartEnd.NumofSeg ++;
			}
			else if (abs(coordend1[0] - coordend2[0]) < 1e-6 && abs(coordend1[1] - coordend2[1]) < 1e-6 && abs(coordend1[2] - coordend2[2]) < 1e-6)
			{
				InsertPartEnd.SegIDs[InsertPartEnd.NumofSeg] = outCellId2;
				InsertPartEnd.SegDir[InsertPartEnd.NumofSeg] = -1;
				InsertPartEnd.NumofSeg ++;
			}
		}
		if (BegAlreadyIn == false)
		{
			if (InsertPartBeg.NumofSeg > 1)
			{
				memcpy(&InsertPart[NumofInsectParts], &InsertPartBeg, sizeof(CInsectPart));
				NumofInsectParts ++;
			}
		}
		if (EndAlreadyIn == false)
		{
			if (InsertPartEnd.NumofSeg > 1)
			{
				memcpy(&InsertPart[NumofInsectParts], &InsertPartEnd, sizeof(CInsectPart));
				NumofInsectParts ++;
			}
		}
	}

	int NumofDeletedVessel = 0;
	int CellIdofDeletedVessel[20];


	int Idx[2];
	CInsectPart MergedInsertPart;  // find very closed bifurcations
/*	for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
	{
		double coord1[3], coord2[3];
		out0Points->GetPoint(pts[0], coord1);
		out0Points->GetPoint(pts[npts-1], coord2);

		if (vtkMath::Distance2BetweenPoints(coord1, coord2) < 5.0*5.0) // too short
		{
			int num = 0;
			for (int i = 0; i < NumofInsectParts; i ++)
			{
				for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
				{
					if(InsertPart[i].SegIDs[j] == outCellId)
					{
						Idx[num] = i;
						num ++;
					}
				}
				if(num >= 2)
				{
					CellIdofDeletedVessel[NumofDeletedVessel] = outCellId;
					NumofDeletedVessel ++;
					break;
				}
			}
			if (num == 2)
			{
				for(int l = 0; l < 3; l ++)
					MergedInsertPart.coord[l] = 0.5 * (InsertPart[Idx[0]].coord[l] + InsertPart[Idx[1]].coord[l]);
				MergedInsertPart.ts = pts[npts/2];
				MergedInsertPart.NumofSeg = InsertPart[Idx[0]].NumofSeg + InsertPart[Idx[1]].NumofSeg - 2;
				int MergedSegIdx = 0;
				for (int k = 0; k < InsertPart[Idx[0]].NumofSeg; k ++)
				{
					if (InsertPart[Idx[0]].SegIDs[k] == outCellId)
						continue;
					MergedInsertPart.SegIDs[MergedSegIdx] = InsertPart[Idx[0]].SegIDs[k];
					MergedInsertPart.SegDir[MergedSegIdx] = InsertPart[Idx[0]].SegDir[k];
					MergedInsertPart.endfaceids[MergedSegIdx] = InsertPart[Idx[0]].endfaceids[k];
					MergedSegIdx ++;
				}
				for (int k = 0; k < InsertPart[Idx[1]].NumofSeg; k ++)
				{
					if (InsertPart[Idx[1]].SegIDs[k] == outCellId)
						continue;
					MergedInsertPart.SegIDs[MergedSegIdx] = InsertPart[Idx[1]].SegIDs[k];
					MergedInsertPart.SegDir[MergedSegIdx] = InsertPart[Idx[1]].SegDir[k];
					MergedInsertPart.endfaceids[MergedSegIdx] = InsertPart[Idx[1]].endfaceids[k];
					MergedSegIdx ++;
				}
				memcpy(&InsertPart[Idx[0]], &MergedInsertPart, sizeof(CInsectPart));
				memcpy(&InsertPart[Idx[1]], &InsertPart[Idx[1] + 1], (NumofInsectParts - Idx[1] - 1) * sizeof(CInsectPart));
				NumofInsectParts --;
			}
		}
	}
*/
//////////////////////////////////////////

	vtkSmartPointer<vtkTriangle> mergeTriangle = vtkSmartPointer<vtkTriangle>::New();
	vtkIdType* StartIDofSegment = new vtkIdType[out0Lines->GetNumberOfCells()];
	vtkIdType* EndIDofSegment = new vtkIdType[out0Lines->GetNumberOfCells()];
	CSegment Segment_circ[20];
	CSegment** Segment_lumn_list = (CSegment**)malloc(NumofInsectParts * sizeof(CSegment*));


	for (int i = 0; i < NumofInsectParts; i ++)
	{
		for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
		{
			for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
			{
				if (outCellId == InsertPart[i].SegIDs[j])
					break;
			}
			
			if (InsertPart[i].SegDir[j] == 1)
			{
				InsertPart[i].endfaceids[j] = 1;
				if (InsertPart[i].endfaceids[j] > npts - 1)
					InsertPart[i].endfaceids[j] = npts - 1;
			}
			else
			{
				InsertPart[i].endfaceids[j] = npts - 2;
				if (InsertPart[i].endfaceids[j] < 0)
					InsertPart[i].endfaceids[j] = 0;
			}
			double coord[3];
			out0Points->GetPoint(pts[InsertPart[i].endfaceids[j]], coord);
			Segment_circ[j].PointList.x = coord[0];
			Segment_circ[j].PointList.y = coord[1];
			Segment_circ[j].PointList.z = coord[2];
		}
	}

	bool isconvexendface[8];

	int LOOPNUM = 0;
	int LOOPNUMbig = 0;

	for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
	{
		StartIDofSegment[outCellId] = 0;
		EndIDofSegment[outCellId] = npts - 1;
	}

	CRingPoint** contour_out_list = (CRingPoint**)malloc(NumofInsectParts * sizeof(CRingPoint*));
	int** PointNuminEachContour_list = (int**)malloc(NumofInsectParts * sizeof(int*));
	int* NumContour_list = (int*)malloc(NumofInsectParts * sizeof(int));

	int NumContour_2 = 0;
	CRingPoint* contour_2 = (CRingPoint*)malloc( 300 * 3 * sizeof(CRingPoint));


	double avertTH = 0.02; 
	double ahoriTH = 0.01;

	for (int i = 0; i < NumofInsectParts; i ++)
	{
		contour_out_list[i] = (CRingPoint*)malloc(MaxNumofContour * MaxNumofPointsInEachContour * sizeof(CRingPoint));
		PointNuminEachContour_list[i] = (int*)malloc(MaxNumofContour * sizeof(int));
		Segment_lumn_list[i] = (CSegment*)malloc(20 * sizeof(CSegment));
		EndFaceisfind[i] = false;
		for(int j = 0; j < InsertPart[i].NumofSeg; j ++)
		{
			isconvexendface[j] = true;
		}
	}

/*******************************************************/
	double* VesselRadius = new double[NumofInsectParts * 20];

	bool* flag_Radius_of_bifur_isfound = new bool[NumofInsectParts];
	double* Radius_of_bifur = new double[NumofInsectParts];

	
	for (int i = 0; i < NumofInsectParts; i ++)
	{ 
		for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
		{
			VesselRadius[20*i+j] = 1.0;
		}

		double vesselradius[20];
		double Radius_of_bifur_initial = 2.0;

FINDBIFURRADIUS:
		for (Radius_of_bifur[i] = Radius_of_bifur_initial; Radius_of_bifur[i] < 8.0; Radius_of_bifur[i] += 0.1)
		{
			flag_Radius_of_bifur_isfound[i] = true;

			for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
			{
				for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
				{
					if (outCellId == InsertPart[i].SegIDs[j])
						break;
				}
				double clcoord[3];
				double distance_choosen = 100.0;
				if (InsertPart[i].SegDir[j] == 1)
				{
					for (int k = 1; k < npts; k ++)
					{
						out0Points->GetPoint(pts[k], clcoord);
						double distance = sqrt(vtkMath::Distance2BetweenPoints(clcoord, InsertPart[i].coord));
						if (abs(distance - Radius_of_bifur[i]) < abs(distance_choosen - Radius_of_bifur[i]))
						{
							distance_choosen = distance;
							InsertPart[i].endfaceids[j] = k;
							Segment_circ[j].PointList.x = clcoord[0];
							Segment_circ[j].PointList.y = clcoord[1];
							Segment_circ[j].PointList.z = clcoord[2];
							Segment_lumn_list[i][j].PointList.x = clcoord[0];
							Segment_lumn_list[i][j].PointList.y = clcoord[1];
							Segment_lumn_list[i][j].PointList.z = clcoord[2];
						}
						if ( abs(distance_choosen - Radius_of_bifur[i]) < 0.1 )
							break;
					}
				}
				else
				{
					for (int k = npts - 2; k >= 0; k --)
					{
						out0Points->GetPoint(pts[k], clcoord);
						double distance = sqrt(vtkMath::Distance2BetweenPoints(clcoord, InsertPart[i].coord));
						if (abs(distance - Radius_of_bifur[i]) < abs(distance_choosen - Radius_of_bifur[i]))
						{
							distance_choosen = distance;
							InsertPart[i].endfaceids[j] = k;
							Segment_circ[j].PointList.x = clcoord[0];
							Segment_circ[j].PointList.y = clcoord[1];
							Segment_circ[j].PointList.z = clcoord[2];
							Segment_lumn_list[i][j].PointList.x = clcoord[0];
							Segment_lumn_list[i][j].PointList.y = clcoord[1];
							Segment_lumn_list[i][j].PointList.z = clcoord[2];							
						}
						if ( abs(distance_choosen - Radius_of_bifur[i]) < 0.1 )
							break;
					}
				}

				out0LumenRadius->GetTuple(pts[InsertPart[i].endfaceids[j]], refineradii);

				vesselradius[j] = 0.0;
				for (int l = 0; l < clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1); l ++)
					vesselradius[j] += refineradii[l];

				vesselradius[j] = vesselradius[j] / (clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1));

				if (vesselradius[j] >= Radius_of_bifur[i])
				{
					flag_Radius_of_bifur_isfound[i] = false;
					break;
				}
			}

			if (flag_Radius_of_bifur_isfound[i] == false)
				continue;

			flag_Radius_of_bifur_isfound[i] = true;
			for (int j1 = 0; j1 < InsertPart[i].NumofSeg - 1; j1 ++)
			{
				double dirj1[3], dirj2[3];
				double anglej1, anglej2;
				dirj1[0] = Segment_circ[j1].PointList.x - InsertPart[i].coord[0];
				dirj1[1] = Segment_circ[j1].PointList.y - InsertPart[i].coord[1];
				dirj1[2] = Segment_circ[j1].PointList.z - InsertPart[i].coord[2];
				vtkMath::Normalize(dirj1);
				anglej1 = asin(vesselradius[j1] / Radius_of_bifur[i]);
				for (int j2 = j1 + 1; j2 < InsertPart[i].NumofSeg; j2 ++)
				{
					dirj2[0] = Segment_circ[j2].PointList.x - InsertPart[i].coord[0];
					dirj2[1] = Segment_circ[j2].PointList.y - InsertPart[i].coord[1];
					dirj2[2] = Segment_circ[j2].PointList.z - InsertPart[i].coord[2];
					vtkMath::Normalize(dirj2);

					anglej2 = asin(vesselradius[j2] / Radius_of_bifur[i]);
					double angle_j1j2 = acos(vtkMath::Dot(dirj1, dirj2));

					if (angle_j1j2 < anglej1 + anglej2 + 2.0*M_PI/180.0)
					{
						flag_Radius_of_bifur_isfound[i] = false;
						break;
					}
				}
				if (flag_Radius_of_bifur_isfound[i] == false)
					break;
			}

			if (flag_Radius_of_bifur_isfound[i] == true)
			{
				for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
				{
					for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
					{
						if (outCellId == InsertPart[i].SegIDs[j])
							break;
					}
					out0Axis1->GetTuple(pts[InsertPart[i].endfaceids[j]], axis1);
					out0Axis2->GetTuple(pts[InsertPart[i].endfaceids[j]], axis2); 
					out0LumenRadius->GetTuple(pts[InsertPart[i].endfaceids[j]], refineradii);

					for (int k = 0; k < clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1); k ++)
					{
						coord[0] = Segment_circ[j].PointList.x + vesselradius[j] * this->RadiusScale* ( cos(k*cirstep) * axis1[0] + sin(k*cirstep) * axis2[0] );
						coord[1] = Segment_circ[j].PointList.y + vesselradius[j] * this->RadiusScale* ( cos(k*cirstep) * axis1[1] + sin(k*cirstep) * axis2[1] );
						coord[2] = Segment_circ[j].PointList.z + vesselradius[j] * this->RadiusScale* ( cos(k*cirstep) * axis1[2] + sin(k*cirstep) * axis2[2] );

						for (int l = 0; l < 3; l ++) dir[l] = coord[l] - InsertPart[i].coord[l];
						vtkMath::Normalize(dir);
						for (int l = 0; l < 3; l ++) coord[l] = InsertPart[i].coord[l] + Radius_of_bifur[i] * dir[l];

						Segment_circ[j].PointList.rx[k] = coord[0];
						Segment_circ[j].PointList.ry[k] = coord[1];
						Segment_circ[j].PointList.rz[k] = coord[2];
					}

					for (int k = 0; k < clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1); k ++)
					{
						coord[0] = Segment_lumn_list[i][j].PointList.x + refineradii[k] * this->RadiusScale* ( cos(k*cirstep) * axis1[0] + sin(k*cirstep) * axis2[0] );
						coord[1] = Segment_lumn_list[i][j].PointList.y + refineradii[k] * this->RadiusScale* ( cos(k*cirstep) * axis1[1] + sin(k*cirstep) * axis2[1] );
						coord[2] = Segment_lumn_list[i][j].PointList.z + refineradii[k] * this->RadiusScale* ( cos(k*cirstep) * axis1[2] + sin(k*cirstep) * axis2[2] );

						Segment_lumn_list[i][j].PointList.rx[k] = coord[0];
						Segment_lumn_list[i][j].PointList.ry[k] = coord[1];
						Segment_lumn_list[i][j].PointList.rz[k] = coord[2];
					}
				}
				for (int j = 0; j < InsertPart[i].NumofSeg; j ++)
				{
					for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
					{
						if (outCellId == InsertPart[i].SegIDs[j])
							break;
					}

					if (InsertPart[i].SegDir[j] == 1)
						StartIDofSegment[outCellId] = InsertPart[i].endfaceids[j];
					else
						EndIDofSegment[outCellId] = InsertPart[i].endfaceids[j];
				}
				break;
			}
		}

		int MergeResult = 1;
		if (flag_Radius_of_bifur_isfound[i])
		{
			MergeResult = MergeAlgorithm(Segment_circ, Segment_lumn_list[i], clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1), InsertPart[i].NumofSeg, InsertPart[i].coord, contour_out_list[i], PointNuminEachContour_list[i], &NumContour_list[i]);
			for (int j = 0; j < InsertPart[i].NumofSeg; j ++) VesselRadius[20*i+j]  = vesselradius[j];
		}

		if (MergeResult == 1) // if merge failed
		{
		//	printf("i = %d, Radius_of_bifur[i] = %.2f, VesselRadius = %.2f\n", i, Radius_of_bifur[i], VesselRadius[20*i+0]);
			Radius_of_bifur_initial = Radius_of_bifur[i] + 0.1;
			if (Radius_of_bifur_initial < 6.0)
				goto FINDBIFURRADIUS;
		}
	}
		

/*******************************************************/

	vtkPoints *out1Points = vtkPoints::New();
	vtkCellArray *out1Strips = vtkCellArray::New();
	vtkDoubleArray *out1Param = vtkDoubleArray::New();
	out1Param->SetName("Param");
	out1Param->SetNumberOfComponents(2);
	vtkDoubleArray *out1Radius = vtkDoubleArray::New();
	out1Radius->SetName("Radius");
	out1Radius->SetNumberOfComponents(1);

	vtkPoints *out2Points = vtkPoints::New();
	vtkCellArray *out2Strips = vtkCellArray::New();
	vtkDoubleArray *out2Param = vtkDoubleArray::New();
	out2Param->SetName("Param");
	out2Param->SetNumberOfComponents(2);
	vtkDoubleArray *out2Radius = vtkDoubleArray::New();
	out2Radius->SetName("Radius");
	out2Radius->SetNumberOfComponents(1);

	vtkSmartPointer<vtkIntArray>isBifurcation = vtkSmartPointer<vtkIntArray>::New();
	isBifurcation->SetName("isBifurcation");
	isBifurcation->SetNumberOfComponents(1);
	vtkDoubleArray* bifurcationPoints = vtkDoubleArray::New();
	bifurcationPoints->SetNumberOfComponents(3);
	bifurcationPoints->SetName("BifurcationPoints");
		
	vtkPoints *out3Points = vtkPoints::New();
	vtkCellArray *out3Strips = vtkCellArray::New();
	vtkDoubleArray *out3Param = vtkDoubleArray::New();
	out3Param->SetName("Param");
	out3Param->SetNumberOfComponents(2);
	vtkDoubleArray *out3Radius = vtkDoubleArray::New();
	out3Radius->SetName("Radius");
	out3Radius->SetNumberOfComponents(1);

	vtkPoints *out4Points = vtkPoints::New();
	vtkCellArray *out4Strips = vtkCellArray::New();
	vtkDoubleArray *out4Param = vtkDoubleArray::New();
	out4Param->SetName("Param");
	out4Param->SetNumberOfComponents(2);
	vtkDoubleArray *out4Radius = vtkDoubleArray::New();
	out4Radius->SetName("Radius");
	out4Radius->SetNumberOfComponents(1);

	output1->GetCellData()->CopyFieldOff("CircumParam");
	output2->GetCellData()->CopyFieldOff("CircumParam");
	output3->GetCellData()->CopyFieldOff("CircumParam");
//	output4->GetCellData()->CopyFieldOff("CircumParam");
	output1->GetCellData()->CopyAllocate(output0->GetCellData());
	output2->GetCellData()->CopyAllocate(output0->GetCellData());
	output3->GetCellData()->CopyAllocate(output0->GetCellData());
	

	for(inCellId=0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts,pts); inCellId++)
	{
//		printf("%d, %d, %d, %d\n", inCellId, pts[StartIDofSegment[inCellId]], pts[EndIDofSegment[inCellId]], pts[EndIDofSegment[inCellId]] - pts[StartIDofSegment[inCellId]]);

		bool thisisadeletedvessel = false;
		for (int i = 0; i < NumofDeletedVessel; i ++)
		{
			if (inCellId == CellIdofDeletedVessel[i])
			{
				thisisadeletedvessel = true;
				break;
			}
		}
		if (thisisadeletedvessel == true)
			continue;


		vtkIdList *idlist1prev = vtkIdList::New();
		vtkIdList *idlist1curr = vtkIdList::New();
		vtkIdList *idlist2prev = vtkIdList::New();
		vtkIdList *idlist2curr = vtkIdList::New();
		vtkIdList *idlist3prev = vtkIdList::New();
		vtkIdList *idlist3curr = vtkIdList::New();

		out0CircumParam->GetTuple(inCellId, circumparam);


//			for(vtkIdType j=0; j<npts; j++)
		for(vtkIdType j = StartIDofSegment[inCellId]; j <= EndIDofSegment[inCellId]; j++)
		{
//				if(j % this->LongitudinalResampleSteps!=0 && j!=0 && j!=npts-1 ) //					continue;
			if(j % this->LongitudinalResampleSteps != 0 && j != StartIDofSegment[inCellId] && j != EndIDofSegment[inCellId])						continue;

			out0Points->GetPoint(pts[j], center);
			out0Dir->GetTuple(pts[j], dir);
			out0Axis1->GetTuple(pts[j], axis1);
			out0Axis2->GetTuple(pts[j], axis2);
			radius = out0Radius->GetValue(pts[j]);
			out0LumenRadius->GetTuple(pts[j], refineradii);
			out0WallThickness->GetTuple(pts[j], refinethickness);
			longiparam = out0LongiParam->GetValue(pts[j]);

			for(int k = 0; k < clLumenRadius->GetNumberOfComponents()*(this->CircumferentialRefineSteps+1); k ++)
			{
				if( k%this->CircumferentialResampleSteps!=0 && k!=0 ) continue;

				for(int l=0; l<3; l++)
					coord[l] = center[l] + radius * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );
				idlist1curr->InsertNextId(out1Points->InsertNextPoint(coord));
				out1Param->InsertNextTuple2(longiparam, circumparam[k]);
				out1Radius->InsertNextValue(radius);


				for(int l=0; l<3; l++)
					coord[l] = center[l] + refineradii[k] * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );

				idlist2curr->InsertNextId(out2Points->InsertNextPoint(coord));

				double temp[3] = {0.0, 0.0, 0.0};
				bifurcationPoints->InsertNextTuple(temp);

				out2Param->InsertNextTuple2(longiparam, circumparam[k]);
				out2Radius->InsertNextValue(refineradii[k]);
				isBifurcation->InsertNextValue(-1);

				for(int l=0; l<3; l++)
					coord[l] = center[l] + (refineradii[k]+refinethickness[k]) * this->RadiusScale * ( cos(k*cirstep)*axis1[l] + sin(k*cirstep)*axis2[l] );

				idlist3curr->InsertNextId(out3Points->InsertNextPoint(coord));
				out3Param->InsertNextTuple2(longiparam, circumparam[k]);
				out3Radius->InsertNextValue(refineradii[k]+refinethickness[k]);
			}
//				if(j!=0)
			if(j !=  StartIDofSegment[inCellId])
			{
				vtkIdType outcellId;
				outcellId = out1Strips->InsertNextCell(2 * (idlist1curr->GetNumberOfIds() + 1));
				output1->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
				for(int k = 0; k < idlist1curr->GetNumberOfIds(); k ++)
				{
					out1Strips->InsertCellPoint(idlist1curr->GetId(k));
					out1Strips->InsertCellPoint(idlist1prev->GetId(k));
				}
				out1Strips->InsertCellPoint(idlist1curr->GetId(0));
				out1Strips->InsertCellPoint(idlist1prev->GetId(0));

			/*	outcellId = out2Strips->InsertNextCell(2 * (idlist2curr->GetNumberOfIds() + 1));
				output2->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
				for(int k=0; k < idlist2curr->GetNumberOfIds(); k++)
				{
					out2Strips->InsertCellPoint(idlist2curr->GetId(k));
					out2Strips->InsertCellPoint(idlist2prev->GetId(k));
				}
				out2Strips->InsertCellPoint(idlist2curr->GetId(0));
				out2Strips->InsertCellPoint(idlist2prev->GetId(0));
			*/

				vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
				for(int k=0; k < idlist2curr->GetNumberOfIds(); k++)
				{
					int k1 = k, k2 = k + 1;
					k2 = k2>=idlist2curr->GetNumberOfIds()?0:k2;
					triangle->GetPointIds()->SetId(0, idlist2curr->GetId(k1));
					triangle->GetPointIds()->SetId(1, idlist2curr->GetId(k2));
					triangle->GetPointIds()->SetId(2, idlist2prev->GetId(k1));
					outcellId = out2Strips->InsertNextCell(triangle);
					output2->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
					triangle->GetPointIds()->SetId(0, idlist2prev->GetId(k1));
					triangle->GetPointIds()->SetId(1, idlist2prev->GetId(k2));
					triangle->GetPointIds()->SetId(2, idlist2curr->GetId(k2));
					outcellId = out2Strips->InsertNextCell(triangle);
					output2->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
				}

				outcellId = out3Strips->InsertNextCell(2 * (idlist3curr->GetNumberOfIds() + 1));
				output3->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
				for(int k = 0; k < idlist3curr->GetNumberOfIds(); k ++)
				{
					out3Strips->InsertCellPoint(idlist3curr->GetId(k));
					out3Strips->InsertCellPoint(idlist3prev->GetId(k));
				}
				out3Strips->InsertCellPoint(idlist3curr->GetId(0));
				out3Strips->InsertCellPoint(idlist3prev->GetId(0));
			}
			idlist1prev->DeepCopy(idlist1curr);
			idlist2prev->DeepCopy(idlist2curr);
			idlist3prev->DeepCopy(idlist3curr);
			idlist1curr->Reset();
			idlist2curr->Reset();
			idlist3curr->Reset();
		}
		idlist1prev->Delete();
		idlist1curr->Delete();
		idlist2prev->Delete();
		idlist2curr->Delete();
		idlist3prev->Delete();
		idlist3curr->Delete();
	}


/*	outcellId = out2Strips->InsertNextCell(2 * (idlist2curr->GetNumberOfIds() + 1));
	output2->GetCellData()->CopyData(output0->GetCellData(),inCellId,outcellId);
	for(int k=0; k < idlist2curr->GetNumberOfIds(); k++)
	{
		out2Strips->InsertCellPoint(idlist2curr->GetId(k));
		out2Strips->InsertCellPoint(idlist2prev->GetId(k));
	}
	out2Strips->InsertCellPoint(idlist2curr->GetId(0));
	out2Strips->InsertCellPoint(idlist2prev->GetId(0));
*/
	
	for (int i = 0; i < NumofInsectParts; i ++)
	{
	//	if (EndFaceisfind[i] == false)
	//		continue;
		if (flag_Radius_of_bifur_isfound[i] == false)
			continue;

		double coordtemp[3];
		vtkIdType existID = 0;
		vtkIdType npts_Before = out2Points->GetNumberOfPoints();

		bool finditinexistOut2points = false;
		for (int j = 0; j < NumContour_list[i]; j ++)
		{
			for (int k = 0; k < PointNuminEachContour_list[i][j]; k ++)
			{
				coord[0] = contour_out_list[i][j * MaxNumofPointsInEachContour + k].x;
				coord[1] = contour_out_list[i][j * MaxNumofPointsInEachContour + k].y;
				coord[2] = contour_out_list[i][j * MaxNumofPointsInEachContour + k].z;
				// find this coord in exist out5Points
				finditinexistOut2points = false;
				if(contour_out_list[i][j * MaxNumofPointsInEachContour + k].l > -1)
				{
					for (vtkIdType kk = 0; kk < out2Points->GetNumberOfPoints(); kk ++)
					{
						out2Points->GetPoint(kk, coordtemp);
						if (abs(coord[0] - coordtemp[0]) < 1e-4
							&& abs(coord[1] - coordtemp[1]) < 1e-4
							&& abs(coord[2] - coordtemp[2]) < 1e-4)
						{
							finditinexistOut2points = true;
							existID = kk;
							break;
						}
					}
				}
				else
				{
					for (vtkIdType kk = npts_Before; kk < out2Points->GetNumberOfPoints(); kk ++)
					{
						out2Points->GetPoint(kk, coordtemp);
						if (abs(coord[0] - coordtemp[0]) < 1e-4
							&& abs(coord[1] - coordtemp[1]) < 1e-4
							&& abs(coord[2] - coordtemp[2]) < 1e-4)
						{
							finditinexistOut2points = true;
							existID = kk;
							break;
						}
					}
				}			
				if (finditinexistOut2points == false)
				{
					mergeTriangle->GetPointIds()->SetId(k, out2Points->InsertNextPoint(coord));
					out2Param->InsertNextTuple2(-1.0, -2.0);
					out2Radius->InsertNextValue(1.0);
					isBifurcation->InsertNextValue(i);

					double temp[3] = {0.0, 0.0, 0.0};
					bifurcationPoints->InsertNextTuple(temp);
				}
				else
					mergeTriangle->GetPointIds()->SetId(k, existID);
			} 

			vtkIdType outcellId = out2Strips->InsertNextCell(mergeTriangle);
			output2->GetCellData()->CopyData(output0->GetCellData(),0,outcellId);  // temp
		}
		bifurcationPoints->SetTuple(i, InsertPart[i].coord);
	}

	output1->SetPoints(out1Points); out1Points->Delete();
	output2->SetPoints(out2Points); 
	output3->SetPoints(out3Points); out3Points->Delete();
	output1->GetPointData()->AddArray(out1Param); out1Param->Delete();
	output2->GetPointData()->AddArray(out2Param); 
	output3->GetPointData()->AddArray(out3Param); out3Param->Delete();				
	output2->GetPointData()->AddArray(isBifurcation);				
	output1->GetPointData()->SetScalars(out1Radius); out1Radius->Delete();
	output2->GetPointData()->SetScalars(out2Radius); 
	output3->GetPointData()->SetScalars(out3Radius); out3Radius->Delete();	
	output1->SetStrips(out1Strips); out1Strips->Delete();
	output2->SetStrips(out2Strips); 
	output3->SetStrips(out3Strips); out3Strips->Delete();

	output2->GetPointData()->AddArray(bifurcationPoints); bifurcationPoints->Delete();


	// refine the convex hull points
	double* vesselcenter = new double[3*20];
	//smoothvtkpolydata_strips(output2, 1000);
	if (1)
	{
//		smoothvtkpolydata(output2, 10, 2);
		out2Points = output2->GetPoints();
		for (int pid = 0; pid < out2Points->GetNumberOfPoints(); pid ++)
		{
			if (isBifurcation->GetValue(pid) == -1)
				continue;

			double ConvexHullCoord[3];
			out2Points->GetPoint(pid, ConvexHullCoord);
			double VesselRadius_thisbifurcation[20];
			for (int j = 0; j < InsertPart[isBifurcation->GetValue(pid)].NumofSeg; j ++)
			{
				for(outCellId = 0, out0Lines->InitTraversal(); out0Lines->GetNextCell(npts, pts); outCellId ++)
				{
					if (outCellId == InsertPart[isBifurcation->GetValue(pid)].SegIDs[j])
						break;
				}				
				double temp[3];
				out0Points->GetPoint(pts[InsertPart[isBifurcation->GetValue(pid)].endfaceids[j]], temp);
				for (int l = 0; l < 3; l ++) vesselcenter[3*j+l] = temp[l];	
				VesselRadius_thisbifurcation[j] = VesselRadius[20*isBifurcation->GetValue(pid)+j];
			}
			double coordnew[3];
			RefineConvexHull(ConvexHullCoord, InsertPart[isBifurcation->GetValue(pid)], Radius_of_bifur[isBifurcation->GetValue(pid)], out0Points, out0Lines, vesselcenter, VesselRadius_thisbifurcation, coordnew);
			
			out2Points->SetPoint(pid, coordnew);
		}

		output2->SetPoints(out2Points); 

		smoothvtkpolydata(output2, 5, 2); // smooth whole mesh
	//	smoothvtkpolydata2(output2, 10);		// just smooth the bifurcation mesh
	}

	printf("ExtendTubeFilter end!\n");




	vtkSmartPointer<vtkLoopSubdivisionFilter> subdivisionFilter2 = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
	vtkSmartPointer<vtkTriangleFilter> triangles2 = vtkSmartPointer<vtkTriangleFilter>::New();
	triangles2->SetInputData(output2);
	triangles2->Update();
	subdivisionFilter2->SetInputConnection(triangles2->GetOutputPort());
	subdivisionFilter2->SetNumberOfSubdivisions(0);
	subdivisionFilter2->Update();
	output2->DeepCopy(subdivisionFilter2->GetOutput());

	



	out2Points->Delete();
	out2Param->Delete();
	out2Radius->Delete();
	out2Strips->Delete();


	delete[] InsertPart;
	delete[] EndFaceisfind;
	delete[] StartIDofSegment;
	delete[] EndIDofSegment;
	delete[] Radius_of_bifur;
	delete[] flag_Radius_of_bifur_isfound;
	delete[] VesselRadius;
	delete[] vesselcenter;

	free(contour_2);

	for (int i = 0; i < NumofInsectParts; i ++)
	{
		free(contour_out_list[i]);
		free(PointNuminEachContour_list[i]);
		free(Segment_lumn_list[i]);
	}
	free(contour_out_list);
	free(PointNuminEachContour_list);
	free(NumContour_list);
	free(Segment_lumn_list);
	}	

	delete[] refineradii_2;
	delete[] refineradii;
	delete[] refinethickness;
	delete[] circumparam;
	out0Points->Delete();
	out0Lines->Delete();
	out0Radius->Delete();
	out0LumenRadius->Delete();
	out0WallThickness->Delete();
	out0Dir->Delete();
	out0Axis1->Delete();
	out0Axis2->Delete();
	out0LongiParam->Delete();
	out0CircumParam->Delete();
	

	return 1;
}

void ExtendTubeFilter::InterpolateRefine(vtkCardinalSpline *spline, double *in, int insize, double *out, int refinesteps)
{
	//if(!spline || !in || !out || insize<2 || refinesteps<0 ) return;

	spline->RemoveAllPoints();
	for(int i=0; i<insize; i++)
	{
		spline->AddPoint(i, in[i]);
	}
	spline->ClosedOn();
	spline->Compute();
	int ct=0;
	for(int i=0; i<insize; i++)
	{
		out[ct++] = spline->Evaluate(i);
		double rs = 1.0/(refinesteps+1);
		for(int j=1; j<=refinesteps; j++)
			out[ct++] = spline->Evaluate(i+j*rs);
	}
}

