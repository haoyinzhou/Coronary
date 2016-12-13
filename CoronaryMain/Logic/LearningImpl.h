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
using namespace cv;

class LearningImpl;
class Learning
{
public:
	Learning();
	~Learning();
	LearningImpl* limpl;
private:
};

class LearningImpl
{
public:
	bool LoadLandmarkClassifiers(int num);
	bool LoadLumenAlongNormalsClassifiers();

	std::vector<CvBoost> lmBoost;
	CvBoost				 lwBoost;

private:
	bool LandMarkLoaded;
	bool LumenWallLoaded;

public:
	LearningImpl();
	~LearningImpl();
};


#endif

#endif //__VTK_WRAP__