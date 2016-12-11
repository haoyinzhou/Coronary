#ifndef __VTK_WRAP__

#ifndef __ExtendTubeFilter_h
#define __ExtendTubeFilter_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include "vtkImageData.h"

#include <iostream>
#include <vector>

class vtkCellArray;
class vtkCellData;
class vtkDataArray;
class vtkFloatArray;
class vtkPointData;
class vtkPoints;
class vtkCardinalSpline;



//class VTKFILTERSCORE_EXPORT ExtendTubeFilter : public vtkPolyDataAlgorithm
class ExtendTubeFilter: public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ExtendTubeFilter,vtkPolyDataAlgorithm);

  static ExtendTubeFilter *New();

  vtkSetClampMacro(LongitudinalRefineSteps, int, 0, 10);
  vtkGetMacro(LongitudinalRefineSteps, int);

  vtkSetClampMacro(CircumferentialRefineSteps, int, 0, 10);
  vtkGetMacro(CircumferentialRefineSteps, int);
  
  vtkSetClampMacro(LongitudinalResampleSteps, int, 1, 10);
  vtkGetMacro(LongitudinalResampleSteps, int);

  vtkSetClampMacro(CircumferentialResampleSteps, int, 1, 10);
  vtkGetMacro(CircumferentialResampleSteps, int);

  vtkSetMacro(RadiusScale, double);
  vtkGetMacro(RadiusScale, double);

  vtkSetMacro(InputImageData, vtkImageData*);
  vtkGetMacro(InputImageData, vtkImageData*);

  virtual void SetUpdateSegment(vtkIdType update);

  void SetWillBuildBifurcationMesh(bool);

  // Description:
  // Turn on/off whether to cap the ends with polygons. Initial value is off.
  vtkSetMacro(Capping,int);
  vtkGetMacro(Capping,int);
  vtkBooleanMacro(Capping,int);

protected:
  ExtendTubeFilter();
  ~ExtendTubeFilter();

  
  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

  // Usual data generation method
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int Capping; //control whether tubes are capped
  int LongitudinalRefineSteps;
  int CircumferentialRefineSteps;
  int LongitudinalResampleSteps;
  int CircumferentialResampleSteps;
  double RadiusScale;
  bool firstSegment;
  vtkImageData *InputImageData;

  vtkIdType UpdateSegment;

  bool WillBuildBifurcationMesh; // = true: will build bifurcation mesh

  // Helper methods
  void InterpolateRefine(vtkCardinalSpline *spline, double *in, int insize, double *out, int refinesteps);

private:
  ExtendTubeFilter(const ExtendTubeFilter&);  // Not implemented.
  void operator=(const ExtendTubeFilter&);  // Not implemented.
  
  vtkPolyData *out0cache;
};

#endif

#endif //__VTK_WRAP__