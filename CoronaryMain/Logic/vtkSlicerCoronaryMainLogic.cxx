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
	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	LumenModel = vtkSmartPointer<vtkPolyData>::New();
	imageData = vtkSmartPointer<vtkImageData>::New();
	memset(NodeOrigin, 0, 3 * sizeof(double));
	for (int l = 0; l < 3; l++) NodeSpaceing[l] = 1.0;
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
void vtkSlicerCoronaryMainLogic::ObserveMRMLScene()
{
	this->Superclass::ObserveMRMLScene();
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

	VolumnNode->GetOrigin(NodeOrigin);
	VolumnNode->GetSpacing(NodeSpaceing);
	NodeOrigin[0] = -NodeOrigin[0];
	NodeOrigin[1] = -NodeOrigin[1];

	imageData_original = VolumnNode->GetImageData();
	imageData->DeepCopy(imageData_original);
	imageData->SetOrigin(NodeOrigin);
	imageData->SetSpacing(NodeSpaceing);

	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
	interpolator->SetInterpolationModeToLinear();
	interpolator->SetOutValue(-3024.0);
	interpolator->Initialize(imageData);

/*	int	 imageDims[3];
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
*/
	DetectLandmarks_core(imageData, learn, landmarks, interpolator);

/*	{
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
*/
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
	vtkSmartPointer<vtkImageData> hessianImage = vtkSmartPointer<vtkImageData>::New();

	GenerateHessianImage(imageData, hessianImage, interpolator);

	std::cout << "GenerateHessianImage done!" << std::endl;
//	SaveVTKImage(imageData, "C:\\work\\Coronary_Slicer\\testdata\\imageData.mha");
//	SaveVTKImage(hessianImage, "C:\\work\\Coronary_Slicer\\testdata\\hessianImage.mha");

	double leftOstium[3], rightOstium[3];
	for (int l = 0; l < 3; l++) leftOstium[l] = landmarks[SmartCoronary::LEFT_CORONARY_OSTIUM][l];
	for (int l = 0; l < 3; l++) rightOstium[l] = landmarks[SmartCoronary::RIGHT_CORONARY_OSTIUM][l];

	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	DetectCenterline_core(imageData, hessianImage, centerlineModel, leftOstium, rightOstium);

/*	vtkSmartPointer<vtkTransform> translation =	vtkSmartPointer<vtkTransform>::New();
	double spaceing[3], origin[3];
	for (int l = 0; l < 3; l++) spaceing[l] = 1 / NodeSpaceing[l];
	for (int l = 0; l < 3; l++) origin[l] = -NodeOrigin[l];
	translation->Scale(spaceing);
	translation->Translate(origin);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetInputData(centerlineModel);
	transformFilter->SetTransform(translation);
	transformFilter->Update();

	std::cout << "NodeOrigins: " << " x: " << NodeOrigin[0] << " y: " << NodeOrigin[1] << " z: " << NodeOrigin[2] << std::endl;
	std::cout << "NodeSpacings: " << " x: " << NodeSpaceing[0] << " y: " << NodeSpaceing[1] << " z: " << NodeSpaceing[2] << std::endl;
	
	centerlineModel = transformFilter->GetOutput();
*/
	//SaveVTKImage(imageData_original, "C:\\work\\Coronary_Slicer\\testdata\\imageData_original.mha");
	//SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel.vtp");

	std::cout << "DetectCenterlinesLogic Done! " << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::DetectLumenLogic(vtkMRMLScalarVolumeNode* VolumnNode, vtkMRMLLinearTransformNode* transformNode)
{
	std::cout << "DetectLumenLogic Begin! " << std::endl;


	std::cout << "DetectLumenLogic Done! " << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::BuildMeshLogic()
{
	std::cout << "BuildMeshLogic Begin! " << std::endl;
	std::cout << "number of cl points: " << centerlineModel->GetPoints()->GetNumberOfPoints() << std::endl;

	//SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel.vtp");
		
	vtkMRMLModelNode* clNode = vtkMRMLModelNode::New();
	vtkMRMLModelDisplayNode* clDisplayNode = vtkMRMLModelDisplayNode::New();
	clNode->SetAndObservePolyData(centerlineModel);
	this->GetMRMLScene()->AddNode(clNode);
	clDisplayNode->SetColor(1, 0, 0);
	clDisplayNode->SetScene(this->GetMRMLScene());
	this->GetMRMLScene()->AddNode(clDisplayNode);
	clNode->SetAndObserveDisplayNodeID(clDisplayNode->GetID());

	vtkSmartPointer<ExtendTubeFilter> centerlineTube = vtkSmartPointer<ExtendTubeFilter>::New();
	centerlineTube->SetInputData(centerlineModel);
	centerlineTube->SetInputImageData(imageData);
	centerlineTube->Update();

	LumenModel = centerlineTube->GetOutput(2);
	std::cout << "number of LumenModel points: " << LumenModel->GetPoints()->GetNumberOfPoints() << std::endl;
	SavePolyData(LumenModel, "C:\\work\\Coronary_Slicer\\testdata\\LumenModel.vtp");

	vtkMRMLModelNode* LumenNode = vtkMRMLModelNode::New();
	vtkMRMLModelDisplayNode* LumenDisplayNode = vtkMRMLModelDisplayNode::New();
	LumenNode->SetAndObservePolyData(LumenModel);
	this->GetMRMLScene()->AddNode(LumenNode);
	LumenDisplayNode->SetColor(0, 0, 1);
	LumenDisplayNode->SetScene(this->GetMRMLScene());
	this->GetMRMLScene()->AddNode(LumenDisplayNode);
	LumenNode->SetAndObserveDisplayNodeID(LumenDisplayNode->GetID());
	
	std::cout << "BuildMeshLogic Done! " << std::endl;
	return true;
}


















