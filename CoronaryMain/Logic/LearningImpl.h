#ifndef __VTK_WRAP__

#ifndef __LearningImpl_h__
#define __LearningImpl_h__

#include <vtkObject.h>
#include "vtkSlicerModuleLogic.h"
#include "vtkSlicerCoronaryMainModuleLogicExport.h"

#include <vector>
#include "opencv2/core/core.hpp"
#include "opencv2/ml/ml.hpp"


using namespace std;

//class VTK_SLICER_CORONARYMAIN_MODULE_LOGIC_EXPORT LearningImpl :
//	public vtkSlicerModuleLogic

class LearningImpl
{
public:
	bool LoadLandmarkClassifiers(int num);
//	std::vector<CvBoost> lmBoost;
	
public:
	LearningImpl();
	~LearningImpl();
};


#endif

#endif //__VTK_WRAP__