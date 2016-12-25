#ifndef __VTK_WRAP__


#ifndef __ImageCurvedReformat_h__
#define __ImageCurvedReformat_h__

#include"vtkImageAlgorithm.h"

class VTK_EXPORT ImageCurvedReformat : public vtkImageAlgorithm
{
public:
   vtkTypeMacro(ImageCurvedReformat,vtkImageAlgorithm);
   static ImageCurvedReformat* New();

   vtkSetMacro(SegmentId, vtkIdType);
   vtkGetMacro(SegmentId, vtkIdType);

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
   ImageCurvedReformat();
   ~ImageCurvedReformat();

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
   ImageCurvedReformat(const ImageCurvedReformat&);  // Not implemented.
   void operator=(const ImageCurvedReformat&);  // Not implemented.

   vtkIdType 	SegmentId;

   double		RadialSpacing;
   int			RadialExtent;

   int			TwistIndex;
   
   int          UpdateImage;
};

#endif //__ImageCurvedReformat_h__


#endif //__VTK_WRAP__
