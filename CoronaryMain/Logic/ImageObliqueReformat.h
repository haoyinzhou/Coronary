#ifndef __VTK_WRAP__


#ifndef __ImageObliqueReformat_h__
#define __ImageObliqueReformat_h__

#include"vtkImageAlgorithm.h"
#include <vector>

class VTK_EXPORT ImageObliqueReformat : public vtkImageAlgorithm
{
public:
   vtkTypeMacro(ImageObliqueReformat,vtkImageAlgorithm);
   static ImageObliqueReformat* New();

   vtkSetMacro(SegmentId, vtkIdType);
   vtkGetMacro(SegmentId, vtkIdType);

   vtkSetMacro(PointId, vtkIdType);
   vtkGetMacro(PointId, vtkIdType);

   vtkSetMacro(RadialSpacing, double);
   vtkGetMacro(RadialSpacing, double);

   vtkSetMacro(RadialExtent, int);
   vtkGetMacro(RadialExtent, int);

   virtual void SetUpdateImage(int update);
   vtkBooleanMacro(UpdateImage,int);
   
   vtkImageData* GetOutput();
   vtkDataObject* GetOutput(int);

protected:
   ImageObliqueReformat();
   ~ImageObliqueReformat();

   virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

   virtual int RequestData(vtkInformation *, vtkInformationVector **, 
	                                         vtkInformationVector *);

   virtual int FillInputPortInformation(int port, vtkInformation *info);
   virtual int FillOutputPortInformation( int, vtkInformation*);
   virtual int RequestInformation(vtkInformation*, vtkInformationVector**, 
	                                                vtkInformationVector*);
private:
   ImageObliqueReformat(const ImageObliqueReformat&);  // Not implemented.
   void operator=(const ImageObliqueReformat&);  // Not implemented.

   vtkIdType		SegmentId;
   vtkIdType		PointId;

   double		RadialSpacing;
   int			RadialExtent;

   int          UpdateImage;
};

#endif //__ImageObliqueReformat_h__



#endif //__VTK_WRAP__
