/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// CoronaryMain Logic includes
#include "vtkSlicerCoronaryMainLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerCoronaryMainLogic);

//----------------------------------------------------------------------------
vtkSlicerCoronaryMainLogic::vtkSlicerCoronaryMainLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerCoronaryMainLogic::~vtkSlicerCoronaryMainLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerCoronaryMainLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

bool vtkSlicerCoronaryMainLogic
::DetectLandmarksLogic(vtkMRMLScalarVolumeNode* VolumnNode)
{
	std::cout << "DetectLandmarksLogic Begin! " << std::endl;

	if (VolumnNode == NULL)
	{
		std::cerr << "VolumnNode is NULL" << std::endl;
		return false;
	}

	imageData = VolumnNode->GetImageData();
	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
	interpolator->SetInterpolationModeToLinear();
	interpolator->SetOutValue(-3024.0);
	interpolator->Initialize(imageData);

	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];
	double imageCenter[3];
	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);
	for (int l = 0; l < 3; l ++) imageCenter[l] = imageOrigins[l] + imageSpacings[l] * imageDims[l] / 2.0;

	std::cout << "Dims: " << " x: " << imageDims[0] << " y: " << imageDims[1] << " z: " << imageDims[2] << std::endl;
	std::cout << "imageSpacings: " << " x: " << imageSpacings[0] << " y: " << imageSpacings[1] << " z: " << imageSpacings[2] << std::endl;
	std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
	std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
	
	DetectLandmarks_core(imageData, learn, landmarks, interpolator);
	for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
	{
		std::cout << landmarks[i][0] << ", " << landmarks[i][1] << ", " << landmarks[i][2] << std::endl;
	}

	std::cout << "DetectLandmarksLogic done!" << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::DetectCenterlinesLogic(vtkMRMLScalarVolumeNode* VolumnNode, vtkMRMLLinearTransformNode* transformNode)
{
	std::cout << "DetectCenterlinesLogic Begin! " << std::endl;

	if (VolumnNode == NULL)
	{
		std::cerr << "VolumnNode is NULL" << std::endl;
		return false;
	}
	if (imageData == NULL)
		imageData = VolumnNode->GetImageData();
	if (hessianImage == NULL)
		hessianImage = vtkSmartPointer<vtkImageData>::New();
	GenerateHessianImage(imageData, hessianImage, interpolator);

	double leftOstium[3], rightOstium[3];
	for (int l = 0; l < 3; l++) leftOstium[l] = landmarks[SmartCoronary::LEFT_CORONARY_OSTIUM][l];
	for (int l = 0; l < 3; l++) rightOstium[l] = landmarks[SmartCoronary::RIGHT_CORONARY_OSTIUM][l];
	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	DetectCenterline_core(imageData, hessianImage, centerlineModel, leftOstium, rightOstium);

	SavePolyData(centerlineModel, "C:\\work\\centerlineModel.vtp");

	std::cout << "DetectCenterlinesLogic Done! " << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::DetectLumenLogic(vtkMRMLScalarVolumeNode* VolumnNode, vtkMRMLLinearTransformNode* transformNode)
{
	std::cout << "DetectLumenLogic Begin! " << std::endl;

	if (VolumnNode == NULL)
	{
		std::cerr << "VolumnNode is NULL" << std::endl;
		return false;
	}

	vtkImageData* imgData = VolumnNode->GetImageData();


	return true;
}

bool vtkSlicerCoronaryMainLogic
::BuildMeshLogic(vtkMRMLScalarVolumeNode* VolumnNode, vtkMRMLLinearTransformNode* transformNode)
{
	std::cout << "BuildLumenLogic Begin! " << std::endl;

	if (VolumnNode == NULL)
	{
		std::cerr << "VolumnNode is NULL" << std::endl;
		return false;
	}

	vtkImageData* imgData = VolumnNode->GetImageData();
	
	return true;
}


















