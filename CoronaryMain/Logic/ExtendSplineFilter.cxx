#include "ExtendSplineFilter.h"

#include "vtkCardinalSpline.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

vtkStandardNewMacro(ExtendSplineFilter);
vtkCxxSetObjectMacro(ExtendSplineFilter,Spline,vtkSpline);

ExtendSplineFilter::ExtendSplineFilter()
{
  this->Subdivide = VTK_SUBDIVIDE_SPECIFIED;
  this->MaximumNumberOfSubdivisions = VTK_INT_MAX;
  this->NumberOfSubdivisions = 100;
  this->Length = 0.1;
  this->GenerateTCoords = VTK_TCOORDS_FROM_NORMALIZED_LENGTH;
  this->TextureLength = 1.0;

  this->Spline = vtkCardinalSpline::New();
  this->TCoordMap = vtkFloatArray::New();
}

ExtendSplineFilter::~ExtendSplineFilter()
{
  if (this->Spline)
    {
    this->Spline->Delete();
    this->Spline = 0;
    }

  if (this->TCoordMap)
    {
    this->TCoordMap->Delete();
    this->TCoordMap = 0;
    }
}

int ExtendSplineFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointData *pd=input->GetPointData();
  vtkPointData *outPD=output->GetPointData();
  vtkCellData *cd=input->GetCellData();
  vtkCellData *outCD=output->GetCellData();
  vtkCellArray *inLines;

  vtkPoints *inPts;
  vtkIdType numLines;
  vtkCellArray *newLines;
  vtkIdType numNewPts, numNewCells;
  vtkPoints *newPts;
  vtkIdType npts=0, *pts=NULL;
  vtkFloatArray *newTCoords=NULL;
  int abort=0;
  vtkIdType inCellId, numGenPts;
  int genTCoords = VTK_TCOORDS_OFF;

  // Check input and initialize
  //
  vtkDebugMacro(<<"Splining polylines");

  if ( !(inPts=input->GetPoints()) || inPts->GetNumberOfPoints() < 1 ||
      !(inLines = input->GetLines()) ||
       (numLines = inLines->GetNumberOfCells()) < 1 )
    {
    return 1;
    }

  if ( this->Spline == NULL )
    {
    vtkWarningMacro(<< "Need to specify a spline!");
    return 1;
    }

  // Create the geometry and topology
  numNewPts = this->NumberOfSubdivisions * numLines;
  newPts = vtkPoints::New();
  newPts->Allocate(numNewPts);
  newLines = vtkCellArray::New();
  newLines->Allocate(newLines->EstimateSize(1,numNewPts));

  // Point data
  if ( (this->GenerateTCoords == VTK_TCOORDS_FROM_SCALARS &&
        pd->GetScalars() != NULL) ||
       (this->GenerateTCoords == VTK_TCOORDS_FROM_LENGTH ||
        this->GenerateTCoords == VTK_TCOORDS_FROM_NORMALIZED_LENGTH) )
    {
    genTCoords = this->GenerateTCoords;
    newTCoords = vtkFloatArray::New();
    newTCoords->SetNumberOfComponents(2);
    newTCoords->Allocate(numNewPts);
    newTCoords->SetName("TCoords");
    outPD->CopyTCoordsOff();
    }
  outPD->InterpolateAllocate(pd,numNewPts);
  this->TCoordMap->Allocate(VTK_CELL_SIZE);

  // Copy cell data
  numNewCells = numLines;
  outCD->CopyNormalsOff();
  outCD->CopyAllocate(cd,numNewCells);

  // Set up the splines
  this->XSpline = this->Spline->NewInstance();
  this->XSpline->DeepCopy(this->Spline);
  this->YSpline = this->Spline->NewInstance();
  this->YSpline->DeepCopy(this->Spline);
  this->ZSpline = this->Spline->NewInstance();
  this->ZSpline->DeepCopy(this->Spline);

  std::map<vtkIdType, vtkIdType> idmap;
  //  Create points along each polyline.
  //
  for (inCellId=0, inLines->InitTraversal();
       inLines->GetNextCell(npts,pts) && !abort; inCellId++)
    {
      this->UpdateProgress(static_cast<double>(inCellId)/numLines);
    abort = this->GetAbortExecute();

    if (npts < 2)
      {
      vtkWarningMacro(<< "Less than two points in line!");
      continue; //skip tubing this polyline
      }

    // Generate the points around the polyline. The strip is not created
    // if the polyline is bad.
    //
    this->TCoordMap->Reset();
	vtkIdList *idlist = vtkIdList::New();
    numGenPts = this->GeneratePoints(idlist, npts, pts, inPts, newPts,
                                     pd, outPD, genTCoords, newTCoords, idmap);
    if ( ! numGenPts )
      {
      //vtkWarningMacro(<< "Could not generate points!");
      continue; //skip splining
      }

    // Generate the polyline
    //
    this->GenerateLine(idlist,numGenPts,inCellId,cd,outCD,newLines);

	idlist->Delete();
    }//for all polylines

  // Update ourselves
  //
  this->TCoordMap->Initialize();

  this->XSpline->Delete();
  this->YSpline->Delete();
  this->ZSpline->Delete();

  output->SetPoints(newPts);
  newPts->Delete();

  output->SetLines(newLines);
  newLines->Delete();

  if ( newTCoords )
    {
    outPD->SetTCoords(newTCoords);
    newTCoords->Delete();
    }

  output->Squeeze();

  return 1;
}

int ExtendSplineFilter::GeneratePoints(vtkIdList *idlist, vtkIdType npts,
                                    vtkIdType *pts, vtkPoints *inPts,
                                    vtkPoints *newPts, vtkPointData *pd,
                                    vtkPointData *outPD, int genTCoords,
                                    vtkFloatArray *newTCoords, 
									std::map<vtkIdType, vtkIdType>& idmap)
{
  vtkIdType i;

  // Initialize the splines
  this->XSpline->RemoveAllPoints();
  this->YSpline->RemoveAllPoints();
  this->ZSpline->RemoveAllPoints();

  // Compute the length of the resulting spline
  double xPrev[3], x[3], length=0.0, len, t, tc, dist;
  inPts->GetPoint(pts[0],xPrev);
  for (i=1; i < npts; i++)
    {
    inPts->GetPoint(pts[i],x);
    len = sqrt(vtkMath::Distance2BetweenPoints(x,xPrev));
    length += len;
    xPrev[0]=x[0]; xPrev[1]=x[1]; xPrev[2]=x[2];
    }
  if ( length <= 0.0 )
    {
    return 0; //failure
    }

  // Now we insert points into the splines with the parametric coordinate
  // based on (polyline) length. We keep track of the parametric coordinates
  // of the points for later point interpolation.
  inPts->GetPoint(pts[0],xPrev);
  for (len=0,i=0; i < npts; i++)
    {
    inPts->GetPoint(pts[i],x);
    dist = sqrt(vtkMath::Distance2BetweenPoints(x,xPrev));
    if (i > 0 && dist == 0)
      {
      continue;
      }
    len += dist;
    t = len/length;
    this->TCoordMap->InsertValue(i,t);

    this->XSpline->AddPoint(t,x[0]);
    this->YSpline->AddPoint(t,x[1]);
    this->ZSpline->AddPoint(t,x[2]);

    xPrev[0]=x[0]; xPrev[1]=x[1]; xPrev[2]=x[2];
    }

  // Compute the number of subdivisions
  vtkIdType numDivs, numNewPts;
  if ( this->Subdivide == VTK_SUBDIVIDE_SPECIFIED )
    {
    numDivs = this->NumberOfSubdivisions;
    }
  else
    {
    numDivs = static_cast<int>(length / this->Length);
    }
  numDivs = ( numDivs < 1 ? 1 : (numDivs > this->MaximumNumberOfSubdivisions ?
                                 this->MaximumNumberOfSubdivisions : numDivs));

  // Now compute the new points
  numNewPts = numDivs + 1;
  vtkIdType idx;
  double s, s0=0.0;
  if ( genTCoords == VTK_TCOORDS_FROM_SCALARS )
    {
    s0=pd->GetScalars()->GetTuple1(pts[0]);
    }
  double tLo = this->TCoordMap->GetValue(0);
  double tHi = this->TCoordMap->GetValue(1);
  for (idx=0, i=0; i < numNewPts; i++)
    {
      t = static_cast<double>(i) / numDivs;
    x[0] = this->XSpline->Evaluate(t);
    x[1] = this->YSpline->Evaluate(t);
    x[2] = this->ZSpline->Evaluate(t);
	if( i==0 )
	{
		if( idmap.find(pts[0]) != idmap.end() ) idlist->InsertNextId(idmap[pts[0]]);
		else
		{
			vtkIdType newId = newPts->InsertNextPoint(x);
			idlist->InsertNextId(newId);
			idmap[pts[0]] = newId;
		}
	}
	else if( i==numNewPts-1 )
	{
		if( idmap.find(pts[npts-1]) != idmap.end() ) idlist->InsertNextId(idmap[pts[npts-1]]);
		else
		{
			vtkIdType newId = newPts->InsertNextPoint(x);
			idlist->InsertNextId(newId);
			idmap[pts[npts-1]] = newId;
		}
	}
	else
	{
		vtkIdType newId = newPts->InsertNextPoint(x);
		idlist->InsertNextId(newId);
	}

    // interpolate point data
    while ( t > tHi && idx < (npts-2))
      {
      idx++;
      tLo = this->TCoordMap->GetValue(idx);
      tHi = this->TCoordMap->GetValue(idx+1);
      }
    tc = (t - tLo) / (tHi - tLo);
	outPD->InterpolateEdge(pd,idlist->GetId(i),pts[idx],pts[idx+1],tc);

    // generate texture coordinates if desired
    if ( genTCoords != VTK_TCOORDS_OFF )
      {
      if ( genTCoords == VTK_TCOORDS_FROM_NORMALIZED_LENGTH )
        {
        tc = t;
        }
      else if ( genTCoords == VTK_TCOORDS_FROM_LENGTH )
        {
        tc = t * length / this->TextureLength;
        }
      else if ( genTCoords == VTK_TCOORDS_FROM_SCALARS )
        {
        s = outPD->GetScalars()->GetTuple1(idlist->GetId(i)); //data just interpolated
        tc = (s - s0) / this->TextureLength;
        }
      newTCoords->InsertTuple2(idlist->GetId(i),tc,0.0);
      } //if generating tcoords
    } //for all new points

  return numNewPts;
}

void ExtendSplineFilter::GenerateLine(vtkIdList *idlist, vtkIdType npts,
                                   vtkIdType inCellId,
                                   vtkCellData *cd, vtkCellData *outCD,
                                   vtkCellArray *newLines)
{
  vtkIdType i, outCellId;

  outCellId = newLines->InsertNextCell(npts);
  outCD->CopyData(cd,inCellId,outCellId);
  for (i=0; i < npts; i++)
    {
		newLines->InsertCellPoint(idlist->GetId(i));
    }
}


const char *ExtendSplineFilter::GetSubdivideAsString()
{
  if ( this->Subdivide == VTK_SUBDIVIDE_SPECIFIED )
    {
    return "Specified by Number of Subdivisions";
    }
  else
    {
    return "Specified by Length";
    }
}

// Return the method of generating the texture coordinates.
const char *ExtendSplineFilter::GetGenerateTCoordsAsString(void)
{
  if ( this->GenerateTCoords == VTK_TCOORDS_OFF )
    {
    return "GenerateTCoordsOff";
    }
  else if ( this->GenerateTCoords == VTK_TCOORDS_FROM_SCALARS )
    {
    return "GenerateTCoordsFromScalar";
    }
  else if ( this->GenerateTCoords == VTK_TCOORDS_FROM_LENGTH )
    {
    return "GenerateTCoordsFromLength";
    }
  else
    {
    return "GenerateTCoordsFromNormalizedLength";
    }
}

void ExtendSplineFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Subdivide: :" << this->GetSubdivideAsString() << "\n";
  os << indent << "Maximum Number of Subdivisions: "
     << this->MaximumNumberOfSubdivisions << "\n";
  os << indent << "Number of Subdivisions: "
     << this->NumberOfSubdivisions << "\n";
  os << indent << "Length: " << this->Length << "\n";
  os << indent << "Spline: " << this->Spline << "\n";
  os << indent << "Generate TCoords: "
     << this->GetGenerateTCoordsAsString() << endl;
  os << indent << "Texture Length: " << this->TextureLength << endl;
}

