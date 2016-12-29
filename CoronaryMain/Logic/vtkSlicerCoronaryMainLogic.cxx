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
	imageData = vtkSmartPointer<vtkImageData>::New();
	imageData_original = vtkSmartPointer<vtkImageData>::New();
	hessianImage = vtkSmartPointer<vtkImageData>::New();
	interpolator = vtkSmartPointer<vtkImageInterpolator>::New();

	WillBuildBifurcationMesh = false;

	memset(landmarks, 0.0, SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS * 3 * sizeof(double));
	memset(NodeOrigin, 0, 3 * sizeof(double));
	for (int l = 0; l < 3; l++) NodeSpaceing[l] = 1.0;

	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	LumenModel = vtkSmartPointer<vtkPolyData>::New();
	centerlineId = vtkSmartPointer<vtkIdFilter>::New();
	
	centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
	LumenModel_display = vtkSmartPointer<vtkPolyData>::New();

	LandmarksRender = vtkSmartPointer<vtkRenderer>::New();
	clRender = vtkSmartPointer<vtkRenderer>::New();
	LumenRender = vtkSmartPointer<vtkRenderer>::New();
	
	addedclnode.clear();
	addedlandmarknode.clear();
	
	cellid_temp = 0;
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
::GetLandMarksCoord(int index, double coord[3])
{
	if (index >= SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS)
		return false;

	for (int l = 0; l < 3; l++)
		coord[l] = landmarks[index][l];

	return true;
}

bool vtkSlicerCoronaryMainLogic
::SetLandMarksCoord(int index, double coord[3])
{
	if (index >= SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS)
		return false;

	for (int l = 0; l < 3; l++)
		this->landmarks[index][l] =	coord[l];

	return true;
}



bool vtkSlicerCoronaryMainLogic
::DetectLandmarksLogic(vtkMRMLScalarVolumeNode* VolumeNode, QProgressBar* progressbar)
{
	std::cout << "DetectLandmarksLogic Begin! " << std::endl;

	if (VolumeNode == NULL)
	{
		std::cerr << "VolumnNode is NULL" << std::endl;
		return false;
	}

	//-------// 
	progressbar->setValue(1);

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
	DetectLandmarks_core(imageData, learn, landmarks, interpolator, progressbar);

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
	
	//-------// 
	progressbar->setValue(100);

	return true;
}

bool vtkSlicerCoronaryMainLogic
::DetectCenterlinesLogic(QProgressBar* progressbar)
{
	std::cout << "DetectCenterlinesLogic Begin! " << std::endl;

	//-------// 
	progressbar->setValue(1);

	hessianImage = vtkSmartPointer<vtkImageData>::New();
	GenerateHessianImage(imageData, hessianImage, interpolator, progressbar);
	std::cout << "GenerateHessianImage done!" << std::endl;

	//-------//
	progressbar->setValue(50);

	double leftOstium[3], rightOstium[3];
	for (int l = 0; l < 3; l++) leftOstium[l] = landmarks[SmartCoronary::LEFT_CORONARY_OSTIUM][l];
	for (int l = 0; l < 3; l++) rightOstium[l] = landmarks[SmartCoronary::RIGHT_CORONARY_OSTIUM][l];

	centerlineModel = vtkSmartPointer<vtkPolyData>::New();
	DetectCenterline_core(imageData, hessianImage, centerlineModel, leftOstium, rightOstium, progressbar);

//	SaveVTKImage(imageData_original, "C:\\work\\Coronary_Slicer\\testdata\\imageData_original.mha");
//	SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel.vtp");

	std::cout << "DetectCenterlinesLogic Done! " << std::endl;

	//-------// 
	progressbar->setValue(100);

	return true;
}

bool vtkSlicerCoronaryMainLogic
::DetectLumenLogic(QProgressBar* progressbar)
{
	std::cout << "DetectLumenLogic Begin! " << std::endl;

	if (imageData == 0)
	{
		std::cerr << "DetectLumenLogic: cannot find image data" << std::endl;
		return false;
	}
	if (centerlineModel->GetNumberOfCells() == 0)
	{
		std::cerr << "cannot find centerline" << std::endl;
		return false;
	}
	std::cout << "DetectLumenLogic has all input data, begin to detect..." << std::endl;


	vtkSmartPointer<vtkIdTypeArray> cidarray = vtkIdTypeArray::SafeDownCast(centerlineModel->GetCellData()->GetArray("SegmentId"));
	for (vtkIdType i = 0; i < cidarray->GetNumberOfTuples(); i++)
	{		
		std::cout << "vessel id = " << i << std::endl;

		if (DetectCenterlineLumenWall_core(centerlineModel, cidarray->GetValue(i), interpolator, learn))
		{
			centerlineModel->Modified();
		}
		progressbar->setValue(100 * (i + 1) / cidarray->GetNumberOfTuples());
	}

	std::cout << "DetectLumenLogic Done! " << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::BuildCenterlinesMeshLogic()
{
	std::cout << "BuildCenterlinesMeshLogic Begin! " << std::endl;
	if (centerlineModel->GetNumberOfCells() == 0)
	{
		std::cerr << "BuildCenterlinesMeshLogic: cannot find centerline model!" << std::endl;
		return false;
	}	

	if (addedclnode.size() != 0)
	{
		for (int i = 0; i < addedclnode.size(); i++)
			this->GetMRMLScene()->RemoveNode(addedclnode.at(i));
		addedclnode.clear();
	}

	//SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel_BuildCenterlinesMeshLogic.vtp");

	centerlineTube = vtkSmartPointer<ExtendTubeFilter>::New();
	centerlineTube->SetWillBuildBifurcationMesh(this->WillBuildBifurcationMesh);
	centerlineTube->SetInputData(centerlineModel);
	centerlineTube->SetInputImageData(imageData);
	centerlineTube->Update();
	
	LumenModel = centerlineTube->GetOutput(2);

	//SavePolyData(centerlineTube->GetOutput(0), "C:\\work\\Coronary_Slicer\\testdata\\centerlineTubeGetOutput(0)_BuildCenterlinesMeshLogic.vtp");
	//SavePolyData(LumenModel, "C:\\work\\Coronary_Slicer\\testdata\\LumenModel_BuildCenterlinesMeshLogic.vtp");

	centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
	LumenModel_display = vtkSmartPointer<vtkPolyData>::New();
	centerlineModel_display->DeepCopy(centerlineModel);
	LumenModel_display->DeepCopy(LumenModel);

	
	clNode = vtkSmartPointer< vtkMRMLModelNode >::New();
	clDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();
	LumenNode = vtkSmartPointer< vtkMRMLModelNode >::New();
	LumenDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();


	vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	transformMatrix->SetElement(0, 0, -1); transformMatrix->SetElement(0, 1, 0);  transformMatrix->SetElement(0, 2, 0);  transformMatrix->SetElement(0, 3, 0);
	transformMatrix->SetElement(1, 0, 0);  transformMatrix->SetElement(1, 1, -1); transformMatrix->SetElement(1, 2, 0);  transformMatrix->SetElement(1, 3, 0);
	transformMatrix->SetElement(2, 0, 0);  transformMatrix->SetElement(2, 1, 0);  transformMatrix->SetElement(2, 2, 1);  transformMatrix->SetElement(2, 3, 0);
	transformMatrix->SetElement(3, 0, 0);  transformMatrix->SetElement(3, 1, 0);  transformMatrix->SetElement(3, 2, 0);  transformMatrix->SetElement(3, 3, 1);

	this->GetMRMLScene()->SaveStateForUndo();
	vtkMRMLNode* thisaddednode;

	thisaddednode = this->GetMRMLScene()->AddNode(clNode);
	addedclnode.push_back(thisaddednode);
	thisaddednode = this->GetMRMLScene()->AddNode(clDisplayNode);
	addedclnode.push_back(thisaddednode);
	clDisplayNode->SetScene(this->GetMRMLScene());
	clNode->SetScene(this->GetMRMLScene());
	clNode->SetName("centerline model");
	clNode->SetAndObserveDisplayNodeID(clDisplayNode->GetID());
	clNode->Modified();
	clDisplayNode->Modified();

	thisaddednode = this->GetMRMLScene()->AddNode(LumenNode);
	addedclnode.push_back(thisaddednode);
	thisaddednode = this->GetMRMLScene()->AddNode(LumenDisplayNode);
	addedclnode.push_back(thisaddednode);
	LumenDisplayNode->SetScene(this->GetMRMLScene());
	LumenNode->SetScene(this->GetMRMLScene());
	LumenNode->SetName("lumen model");
	LumenNode->SetAndObserveDisplayNodeID(LumenDisplayNode->GetID());
	LumenNode->Modified();
	LumenDisplayNode->Modified();

	clDisplayNode->SetColor(1, 0, 0);
	clDisplayNode->SelectableOn();
	clNode->SetAndObservePolyData(centerlineModel_display);
	clNode->ApplyTransformMatrix(transformMatrix);

	LumenDisplayNode->SetColor(0, 0, 1);
	LumenDisplayNode->SetOpacity(0.5);
	LumenDisplayNode->SetVisibility(1);
	LumenDisplayNode->SetRepresentation(vtkMRMLModelDisplayNode::WireframeRepresentation);
	LumenDisplayNode->SelectableOff();
	LumenNode->SetAndObservePolyData(LumenModel_display);
	LumenNode->ApplyTransformMatrix(transformMatrix);
	
	std::cout << "BuildCenterlinesMeshLogic Done! " << std::endl;
	return true;
}

bool vtkSlicerCoronaryMainLogic
::BuildLandmarksMeshLogic()
{
	std::cout << "BuildLandmarksMeshLogic Begin! " << std::endl;

	if (addedlandmarknode.size() != 0)
	{
		for (int i = 0; i < addedlandmarknode.size(); i++)
			this->GetMRMLScene()->RemoveNode(addedlandmarknode.at(i));
		addedlandmarknode.clear();
	}

	// show landmarks
	vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
	for (int i = SmartCoronary::LEFT_CORONARY_OSTIUM; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
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

	vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	transformMatrix->SetElement(0, 0, -1); transformMatrix->SetElement(0, 1, 0);  transformMatrix->SetElement(0, 2, 0);  transformMatrix->SetElement(0, 3, 0);
	transformMatrix->SetElement(1, 0, 0);  transformMatrix->SetElement(1, 1, -1); transformMatrix->SetElement(1, 2, 0);  transformMatrix->SetElement(1, 3, 0);
	transformMatrix->SetElement(2, 0, 0);  transformMatrix->SetElement(2, 1, 0);  transformMatrix->SetElement(2, 2, 1);  transformMatrix->SetElement(2, 3, 0);
	transformMatrix->SetElement(3, 0, 0);  transformMatrix->SetElement(3, 1, 0);  transformMatrix->SetElement(3, 2, 0);  transformMatrix->SetElement(3, 3, 1);

	this->GetMRMLScene()->SaveStateForUndo();
	vtkMRMLNode* thisaddednode;

	LandmarkNode = vtkSmartPointer< vtkMRMLModelNode >::New();
	LandmarkDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();
	thisaddednode = this->GetMRMLScene()->AddNode(LandmarkNode);
	addedlandmarknode.push_back(thisaddednode);
	thisaddednode = this->GetMRMLScene()->AddNode(LandmarkDisplayNode);
	addedlandmarknode.push_back(thisaddednode);

	LandmarkDisplayNode->SetScene(this->GetMRMLScene());
	LandmarkNode->SetScene(this->GetMRMLScene());
	LandmarkNode->SetName("landmarks");
	LandmarkNode->SetAndObserveDisplayNodeID(LandmarkDisplayNode->GetID());
	LandmarkNode->Modified();
	LandmarkDisplayNode->Modified();

	LandmarkDisplayNode->SetColor(1, 0, 0);
	LandmarkNode->SetAndObservePolyData(cleanFilter->GetOutput());
	LandmarkNode->ApplyTransformMatrix(transformMatrix);
	
	std::cout << "BuildCenterlinesMeshLogic Done! " << std::endl;

	return true;
}


class CMyvtkCommand : public vtkCommand
{
public:
	vtkTypeMacro(CMyvtkCommand, vtkCommand);

	static CMyvtkCommand *New()
	{
		return new CMyvtkCommand;
	}

	void Execute(vtkObject *vtkNotUsed(caller), unsigned long vtkNotUsed(eventId),
		void *vtkNotUsed(callData))
	{
		double ras[3];
		this->crosshairNode->GetCursorPositionRAS(ras);
		std::cout << ras[0] << ", " << ras[1] << ", " << ras[2] << std::endl;
	}

public:
	vtkMRMLCrosshairNode* crosshairNode;
};

class CenterlineMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static CenterlineMouseInteractorStyle* New();
	vtkTypeMacro(CenterlineMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

	CenterlineMouseInteractorStyle()
	{

	}

	virtual void OnRightButtonDown()
	{
		int pickPosition[2];
		this->GetInteractor()->GetEventPosition(pickPosition);
		std::cout << pickPosition[0] << ", " << pickPosition[1] << std::endl;
	}

public:
	
};
vtkStandardNewMacro(CenterlineMouseInteractorStyle);




bool vtkSlicerCoronaryMainLogic
::TestLogic()
{
	std::cout << "TestLogic begin" << std::endl;
//	QVesselEditingWidget* temp; // = new QVesselEditingWidget;
//	temp->resize(600, 500);
//	temp->show();

	//emit visibilitychanged(true);

	/*	QVesselEditingWidget *VesselEditingWidget = new QVesselEditingWidget;

	VesselEditingWidget->mainlogic = this;
	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
	qMRMLThreeDView* threeDView = layoutManager->threeDWidget(0)->threeDView();
	VesselEditingWidget->SlicerThreeDWidget = threeDView->VTKWidget();
	VesselEditingWidget->resize(600, 1200);

	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(5.0);
	sphereSource->Update();

	vtkSmartPointer<vtkPolyDataMapper> VesselEditingMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	VesselEditingMapper->SetInputConnection(sphereSource->GetOutputPort());

	vtkSmartPointer<vtkActor> VesselEditingActor = vtkSmartPointer<vtkActor>::New();
	VesselEditingActor->SetMapper(VesselEditingMapper);

	vtkSmartPointer<vtkRenderer> VesselEditingRenderer = vtkSmartPointer<vtkRenderer>::New();
	VesselEditingRenderer->SetBackground(0, 0, 0); // Background color black
	VesselEditingRenderer->AddActor(VesselEditingActor);

	VesselEditingWidget->GetRenderWindow()->AddRenderer(VesselEditingRenderer);
	
	VesselEditingWidget->show();
*/

	

/*	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
	qMRMLThreeDView* threeDView = layoutManager->threeDWidget(0)->threeDView();
	threeDView->cornerAnnotation()->SetText(0, "Haoyin Zhou Test Button!");
	threeDView->cornerAnnotation()->GetTextProperty()->SetColor(1, 0, 0);
	threeDView->forceRender();

	vtkRenderWindowInteractor* RenderWindowInteractorthreeD = threeDView->VTKWidget()->GetInteractor();
	vtkSmartPointer<CenterlineMouseInteractorStyle> style = vtkSmartPointer<CenterlineMouseInteractorStyle>::New();
	RenderWindowInteractorthreeD->SetInteractorStyle(style);
	RenderWindowInteractorthreeD->Initialize();
	RenderWindowInteractorthreeD->Start();
*/
	
/*	vtkSmartPointer<vtkSphereSource> sphereSource =	vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(5.0);
	sphereSource->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =	vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(sphereSource->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =	vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1, 1, 1); // Background color white
	renderer->AddActor(actor);

	vtkSmartPointer<vtkRenderWindow> renderWindow =	vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(300, 1000); //(width, height)
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	renderWindowInteractor->SetInteractorStyle(style);

	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();
*/


//	vtkMRMLCrosshairNode* crosshairNode = NULL;
//	crosshairNode = vtkMRMLCrosshairNode::SafeDownCast(this->GetMRMLScene()->GetNthNodeByClass(0, "vtkMRMLCrosshairNode"));
//	vtkSmartPointer<CMyvtkCommand> tempcommand = vtkSmartPointer<CMyvtkCommand>::New();
//	tempcommand->crosshairNode = crosshairNode;
//	crosshairNode->AddObserver(vtkMRMLCrosshairNode::CursorPositionModifiedEvent, tempcommand);
	

/*	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
qMRMLThreeDView* threeDView = layoutManager->threeDWidget(0)->threeDView();
threeDView->cornerAnnotation()->SetText(0, "Ready for 3D vessel pick..");
threeDView->cornerAnnotation()->GetTextProperty()->SetColor(1, 0, 0);
threeDView->forceRender();

vtkRenderWindowInteractor* RenderWindowInteractorthreeD = threeDView->VTKWidget()->GetInteractor();
vtkSmartPointer<CenterlineMouseInteractorStyle> style = vtkSmartPointer<CenterlineMouseInteractorStyle>::New();
RenderWindowInteractorthreeD->SetInteractorStyle(style);
RenderWindowInteractorthreeD->Initialize();
RenderWindowInteractorthreeD->Start();
*/
	

//	vtkSmartPointer<vtkCollection> nodecollection = 	this->GetMRMLScene()->GetNodes();
//	for (int i = 0; i < nodecollection->GetNumberOfItems(); i++)
//	{
//		vtkSmartPointer<vtkObject> nodeobject =	nodecollection->GetItemAsObject(i);
//		std::cout << nodeobject->GetClassName() << std::endl;
//	}


	std::cout << "TestLogic end" << std::endl;
	return true;
}



















