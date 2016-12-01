#include "LearningImpl.h"
#include <sstream>

using namespace cv;

Learning::Learning()
{
	limpl = new LearningImpl;
}
Learning::~Learning()
{
	delete limpl;
}


LearningImpl::LearningImpl()
{
}

LearningImpl::~LearningImpl()
{
}

bool LearningImpl::LoadLandmarkClassifiers(int num)
{
	std::ostringstream strstm;
	lmBoost.resize(num);
//	std::cout << "Loading the classifiers ..." << std::endl;
	cv::FileStorage fs("C:\\work\\classifiers\\lvcorlmclassifier.yml", cv::FileStorage::READ);
		
	for (int id = 0; id < num; id++)
	{
		strstm.str("");
		strstm << "Classifier_" << id;
		lmBoost[id].read(*fs, *fs[strstm.str().c_str()]);
		if (!lmBoost[id].get_weak_predictors())
		{
			std::cerr << "Could not read " << strstm.str() << std::endl;
			return false;
		}
	}

	fs.release();
	
	std::cout << "classifiers loading done" << std::endl;

	return true;
}
