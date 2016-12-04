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
	memset(landmarks, 0.0, SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS * 3 * sizeof(double));
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

	double NodeOrigin[3];
	double NodeSpaceing[3];
	VolumnNode->GetOrigin(NodeOrigin);
	VolumnNode->GetSpacing(NodeSpaceing);

	imageData = VolumnNode->GetImageData();
	//	imageData->SetOrigin(-NodeOrigin[0], -NodeOrigin[1], NodeOrigin[2]);
	//	imageData->SetSpacing(NodeSpaceing);

	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
	interpolator->SetInterpolationModeToLinear();
	interpolator->SetOutValue(-3024.0);
	interpolator->Initialize(imageData);

	int	 imageDims[3];
	double imageOrigins[3];
	double imageSpacings[3];
	double imageCenter[3];
	double imageBounds[6];
	imageData->GetDimensions(imageDims);
	imageData->GetOrigin(imageOrigins);
	imageData->GetSpacing(imageSpacings);
	imageData->GetBounds(imageBounds);
	for (int l = 0; l < 3; l++) imageCenter[l] = imageOrigins[l] + imageSpacings[l] * imageDims[l] / 2.0;
	std::cout << "Dims: " << " x: " << imageDims[0] << " y: " << imageDims[1] << " z: " << imageDims[2] << std::endl;
	std::cout << "imageSpacings: " << " x: " << imageSpacings[0] << " y: " << imageSpacings[1] << " z: " << imageSpacings[2] << std::endl;
	std::cout << "imageOrigins: " << " x: " << imageOrigins[0] << " y: " << imageOrigins[1] << " z: " << imageOrigins[2] << std::endl;
	std::cout << "imageBounds: " << imageBounds[0] << ", " << imageBounds[1] << ", " << imageBounds[2] << ", " << imageBounds[3] << ", " << imageBounds[4] << ", " << imageBounds[5] << std::endl;
	std::cout << "NodeOrigins: " << " x: " << NodeOrigin[0] << " y: " << NodeOrigin[1] << " z: " << NodeOrigin[2] << std::endl;
	std::cout << "NodeSpacings: " << " x: " << NodeSpaceing[0] << " y: " << NodeSpaceing[1] << " z: " << NodeSpaceing[2] << std::endl;
	std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
	std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;

	DetectLandmarks_core(imageData, learn, landmarks, interpolator);
	//	imageData->SetOrigin(-NodeOrigin[0], -NodeOrigin[1], NodeOrigin[2]);
	//	imageData->SetSpacing(NodeSpaceing);

	SaveVTKImage(imageData, "C:\\work\\Coronary_Slicer\\testdata\\imageData.mha");

	{
		for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
		{
			//	for (int l = 0; l < 3; l++) landmarks[i][l] = landmarks[i][l] * imageSpacings[l] + imageOrigins[l];
			std::cout << "landmarks: " << landmarks[i][0] << ", " << landmarks[i][1] << ", " << landmarks[i][2] << std::endl;
		}

		vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
		for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
		{
			vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->SetCenter(landmarks[i][0], landmarks[i][1], landmarks[i][2]);
			sphereSource->SetRadius(4.0);
			sphereSource->Update();
			appendFilter->AddInputData(sphereSource->GetOutput());
		}
		appendFilter->Update();
		vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
		cleanFilter->Update();
		SavePolyData(cleanFilter->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\landmarks.vtp");
	}

	{
		double landmarks_forsave[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS][3];
		for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
		{
			landmarks_forsave[i][0] = landmarks[i][0] * NodeSpaceing[0] - NodeOrigin[0];
			landmarks_forsave[i][1] = landmarks[i][1] * NodeSpaceing[1] - NodeOrigin[1];
			landmarks_forsave[i][2] = landmarks[i][2] * NodeSpaceing[2] + NodeOrigin[2];
			std::cout << "landmarks_forsave: " << landmarks_forsave[i][0] << ", " << landmarks_forsave[i][1] << ", " << landmarks_forsave[i][2] << std::endl;
		}
		vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
		for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
		{
			vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->SetCenter(landmarks_forsave[i][0], landmarks_forsave[i][1], landmarks_forsave[i][2]);
			sphereSource->SetRadius(4.0);
			sphereSource->Update();
			appendFilter->AddInputData(sphereSource->GetOutput());
		}
		appendFilter->Update();
		vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
		cleanFilter->Update();
		SavePolyData(cleanFilter->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\landmarks_forsave.vtp");
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

	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
	interpolator->SetInterpolationModeToLinear();
	interpolator->SetOutValue(-3024.0);
	interpolator->Initialize(imageData);
	GenerateHessianImage(imageData, hessianImage, interpolator);

	SaveVTKImage(imageData, "C:\\work\\Coronary_Slicer\\testdata\\imageData.mha");
	SaveVTKImage(hessianImage, "C:\\work\\Coronary_Slicer\\testdata\\hessianImage.mha");
	std::cout << "GenerateHessianImage done!" << std::endl;

	double leftOstium[3], rightOstium[3];
	for (int l = 0; l < 3; l++) leftOstium[l] = landmarks[SmartCoronary::LEFT_CORONARY_OSTIUM][l];
	for (int l = 0; l < 3; l++) rightOstium[l] = landmarks[SmartCoronary::RIGHT_CORONARY_OSTIUM][l];
	std::cout << "leftOstium: " << leftOstium[0] << ", " << leftOstium[1] << " , " << leftOstium[2] << std::endl;
	std::cout << "rightOstium: " << rightOstium[0] << ", " << rightOstium[1] << " , " << rightOstium[2] << std::endl;
	
/*	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	DetectCenterline_core(imageData, hessianImage, centerlineModel, leftOstium, rightOstium);

	SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel.vtp");
*/
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


















