#ifndef __VTK_WRAP__


#ifndef __ImageStretchCurvedReformat__
#define __ImageStretchCurvedReformat__

#include"vtkImageAlgorithm.h"

class VTK_EXPORT ImageStretchCurvedReformat : public vtkImageAlgorithm
{
public:
   vtkTypeMacro(ImageStretchCurvedReformat,vtkImageAlgorithm);
   static ImageStretchCurvedReformat* New();

   vtkSetMacro(SegmentId, vtkIdType*);
   vtkGetMacro(SegmentId, vtkIdType*);

   vtkSetMacro(TwistIndex, int);
   vtkGetMacro(TwistIndex, int);

   vtkSetMacro(RadialSpacing, double);
   vtkGetMacro(RadialSpacing, double);

   vtkSetMacro(RadialExtent, int);
   vtkGetMacro(RadialExtent, int);

   virtual void SetUpdateImage(int update);
   vtkBooleanMacro(UpdateImage,int);

   vtkImageData* GetOutput();
   vtkDataObject* GetOutput(int);

protected:
   ImageStretchCurvedReformat();
   ~ImageStretchCurvedReformat();

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
   ImageStretchCurvedReformat(const ImageStretchCurvedReformat&);  // Not implemented.
   void operator=(const ImageStretchCurvedReformat&);  // Not implemented.

   vtkIdType*		SegmentId;

   double		RadialSpacing;
   int			RadialExtent;

   int			TwistIndex;
   
   int          UpdateImage;
};

#endif //__ImageCurvedReformat_h__


#endif //__VTK_WRAP__