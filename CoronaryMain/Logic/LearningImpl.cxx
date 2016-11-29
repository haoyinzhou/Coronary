#include "LearningImpl.h"
#include <sstream>

using namespace cv;

LearningImpl::LearningImpl()
{
	std::cout << "LearningImpl()" << std::endl;
}

LearningImpl::~LearningImpl()
{
	std::cout << "~LearningImpl()" << std::endl;
}

bool LearningImpl::LoadLandmarkClassifiers(int num)
{
//	std::ostringstream strstm;
//	lmBoost.resize(num);
	std::cout << "Loading the classifiers ..." << std::endl;
//	cv::FileStorage fs("C:\\work\\classifiers\\lvcorlmclassifier.yml", cv::FileStorage::READ);

	FileStorage fs;
	
//	fs.open("C:\\work\\classifiers\\lvcorlmclassifier.yml", cv::FileStorage::READ);
/*	for (int id = 0; id < num; id++)
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
*/
//	fs.release();
	
	std::cout << "done" << std::endl;

	return true;
}
