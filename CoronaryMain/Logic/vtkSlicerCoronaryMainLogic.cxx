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
	addedselectedclnode.clear();

	addedctrlobservertag.clear();
	addedvesselpickobservertag.clear();

	cellid_temp = 0;
}

//----------------------------------------------------------------------------
vtkSlicerCoronaryMainLogic::~vtkSlicerCoronaryMainLogic()
{
/*	if (imageData)
		imageData->Delete();
	if (imageData_original)
		imageData_original->Delete();
	if (hessianImage)
		hessianImage->Delete();
	if (interpolator)
		interpolator->Delete();
	if (centerlineModel)
		centerlineModel->Delete();
	if (LumenModel)
		LumenModel->Delete();
	if (centerlineId)
		centerlineId->Delete();
*/
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

	//SavePolyData(centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\centerlineModel.vtp");

	centerlineTube = vtkSmartPointer<ExtendTubeFilter>::New();
	centerlineTube->SetWillBuildBifurcationMesh(this->WillBuildBifurcationMesh);
	centerlineTube->SetInputData(centerlineModel);
	centerlineTube->SetInputImageData(imageData);
	centerlineTube->Update();
	LumenModel = centerlineTube->GetOutput(2);

//	SavePolyData(LumenModel, "C:\\work\\Coronary_Slicer\\testdata\\LumenModel.vtp");

	centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
	LumenModel_display = vtkSmartPointer<vtkPolyData>::New();
	centerlineModel_display->DeepCopy(centerlineModel);
	LumenModel_display->DeepCopy(LumenModel);

	//
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

bool vtkSlicerCoronaryMainLogic
::RemoveAllSelectedVesselThreeD()
{
	if (addedselectedclnode.size() != 0)
	{
		for (int i = 0; i < addedselectedclnode.size(); i++)
			this->GetMRMLScene()->RemoveNode(addedselectedclnode.at(i));
		addedselectedclnode.clear();
	}

	return true;
}

bool vtkSlicerCoronaryMainLogic
::ShowSelectedVesselThreeD(vtkIdType cellid)
{
	std::cout << "ShowSelectedVesselThreeD begin!" << std::endl;

	vtkSmartPointer<vtkIdTypeArray> centerlineSelectId = vtkSmartPointer<vtkIdTypeArray>::New();
	centerlineSelectId->SetName("SegmentId");
	centerlineSelectId->SetNumberOfValues(1);
	centerlineSelectId->SetValue(0, cellid);

	if (cellid >= centerlineModel->GetNumberOfCells())
		cellid = 0;

	RemoveAllSelectedVesselThreeD();

	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::VALUES);
	selectionNode->SetSelectionList(centerlineSelectId);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	//	extractSelection->SetInputConnection(0, centerlineTube->GetOutputPort(2));
	extractSelection->SetInputData(0, LumenModel_display);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();

	vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());

	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	geometryFilter->SetInputData(selected);
	geometryFilter->Update();

	vtkSmartPointer<vtkPolyData> selectpolydata = geometryFilter->GetOutput();


	SelectedClNode = vtkSmartPointer< vtkMRMLModelNode >::New();
	SelectedClDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();

	vtkMRMLNode* thisaddednode;

	thisaddednode = this->GetMRMLScene()->AddNode(SelectedClNode);
	addedselectedclnode.push_back(thisaddednode);
	thisaddednode = this->GetMRMLScene()->AddNode(SelectedClDisplayNode);
	addedselectedclnode.push_back(thisaddednode);
	SelectedClDisplayNode->SetScene(this->GetMRMLScene());
	SelectedClNode->SetScene(this->GetMRMLScene());
	SelectedClNode->SetName("selected vessel");
	SelectedClNode->SetAndObserveDisplayNodeID(SelectedClDisplayNode->GetID());
	SelectedClNode->Modified();
	SelectedClDisplayNode->Modified();

	SelectedClDisplayNode->SetColor(1, 1, 0);
	SelectedClDisplayNode->SetPointSize(5);
	SelectedClDisplayNode->SetLineWidth(5);
	SelectedClDisplayNode->SetVisibility(1);
	SelectedClDisplayNode->LightingOn();
	SelectedClDisplayNode->SetRepresentation(vtkMRMLModelDisplayNode::WireframeRepresentation);
	SelectedClNode->SetAndObservePolyData(selectpolydata);
	
	std::cout << "ShowSelectedVesselThreeD end!" << std::endl;

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


// mouse pick vessel callback

class CVesselPickCallBack : public vtkCommand
{
public:
	static CVesselPickCallBack *New() { return new CVesselPickCallBack; }

	CVesselPickCallBack()
	{
		LastSelectedID = -1;
	}

	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		std::cout << "Mouse click!" << std::endl;
		if (Slicer3DRenderWindowInteractor == NULL)
			return;
		if (clmodel->GetNumberOfCells() == 0)
			return;

		int pickpixel[2];
		Slicer3DRenderWindowInteractor->GetEventPosition(pickpixel);
		//std::cout << "pickpixel = " << pickpixel[0] << ", " << pickpixel[1] << std::endl;

		vtkSmartPointer< vtkCellPicker > picker = vtkCellPicker::SafeDownCast(Slicer3DRenderWindowInteractor->GetPicker());
		picker->Pick(pickpixel[0], pickpixel[1], 0, Slicer3DRender);

		vtkIdType pickid = picker->GetCellId();
		if (pickid == -1)
			return;
	
		vtkSmartPointer<vtkPolyData> pickedpoly = vtkPolyData::SafeDownCast(picker->GetDataSet());
		
		if (pickedpoly->GetNumberOfCells() == lumenmodel->GetNumberOfCells())
		{
			vtkSmartPointer<vtkIdTypeArray> segmentidarray = vtkIdTypeArray::SafeDownCast(lumenmodel->GetCellData()->GetArray("SegmentId"));
			pickid = segmentidarray->GetValue(pickid);
		}

		if (pickid == LastSelectedID)
			return;
		LastSelectedID = pickid;

		std::cout << "pickid = " << pickid << std::endl;

		logic->ShowSelectedVesselThreeD(pickid);
		
		// pop the vessel editing window
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

		vtkSmartPointer<vtkRenderWindow> VesselEditingRenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
		int slicerthreeDwindowsize[2];
		Slicer3DRenderWindowInteractor->GetSize(slicerthreeDwindowsize);
		VesselEditingRenderWindow->SetSize(0.2 * slicerthreeDwindowsize[0], 1.2 * slicerthreeDwindowsize[1]); //(width, height)
		VesselEditingRenderWindow->AddRenderer(VesselEditingRenderer);

		vtkSmartPointer<vtkRenderWindowInteractor> VesselEditingRenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		VesselEditingRenderWindowInteractor->SetRenderWindow(VesselEditingRenderWindow);

		//vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
		//renderWindowInteractor->SetInteractorStyle(style);
		
	//	VesselEditingRenderWindowInteractor->Start(); 
	

	}

public:
	vtkRenderWindowInteractor* Slicer3DRenderWindowInteractor;
	QVTKWidget* Slicer3DWidget;
	vtkRenderer* Slicer3DRender;
	vtkSlicerCoronaryMainLogic* logic;

	vtkPolyData* clmodel;
	vtkPolyData* lumenmodel;

private:
	vtkIdType LastSelectedID;

};


// Ctrl Key Event
class vtkCtrlKeyPressedInteractionCallback : public vtkCommand
{
public:
	static vtkCtrlKeyPressedInteractionCallback *New()
	{
		return new vtkCtrlKeyPressedInteractionCallback;
	}

	vtkCtrlKeyPressedInteractionCallback()
	{

	}

	~vtkCtrlKeyPressedInteractionCallback()
	{

	}

	void SetInteractor(vtkRenderWindowInteractor *iren)
	{
		this->Iren = iren;
	}

	virtual void Execute(vtkObject *caller, unsigned long ev, void *)
	{
		if (clmodel->GetNumberOfCells() == 0)
			return;
		if (vtkStdString(Iren->GetKeySym()) == "Control_L")
		{
			std::cout << "Control_L is pressed!" << std::endl;
			Iren->SetPicker(VesselPicker);
			addedvesselpickobservertag->push_back(Iren->AddObserver(vtkCommand::LeftButtonPressEvent, VesselPickCallBack, 10.0f));
		}
	}

private:
	vtkRenderWindowInteractor* Iren;
public:
	vtkCellPicker* VesselPicker;
	CVesselPickCallBack* VesselPickCallBack;

	vector<unsigned long>* addedvesselpickobservertag;
	vtkPolyData* clmodel;
};

class vtkCtrlKeyReleasedInteractionCallback : public vtkCommand
{
public:
	static vtkCtrlKeyReleasedInteractionCallback *New()
	{
		return new vtkCtrlKeyReleasedInteractionCallback;
	}

	vtkCtrlKeyReleasedInteractionCallback()
	{

	}

	~vtkCtrlKeyReleasedInteractionCallback()
	{

	}

	void SetInteractor(vtkRenderWindowInteractor *iren)
	{
		this->Iren = iren;
	}

	virtual void Execute(vtkObject *caller, unsigned long ev, void *)
	{
		if (vtkStdString(Iren->GetKeySym()) == "Control_L")
		{
			std::cout << "Control_L is released!" << std::endl;
			for (int i = 0; i < addedvesselpickobservertag->size(); i++)
				Iren->RemoveObserver(addedvesselpickobservertag->at(i));
			addedvesselpickobservertag->clear();
		}
	}

private:
	vtkRenderWindowInteractor*	Iren;

public:
	vector<unsigned long>* addedvesselpickobservertag;

};


bool vtkSlicerCoronaryMainLogic
::SetupKeyMouseObserver()
{
	if (centerlineModel->GetNumberOfCells() == 0)
		return false;

	RemoveAllSelectedVesselThreeD();

	// set ctrl observer
	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
	QVTKWidget* threeDView = layoutManager->threeDWidget(0)->threeDView()->VTKWidget();
	vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractorthreeD = threeDView->GetInteractor();
	vtkSmartPointer<vtkRendererCollection> rendercollection = threeDView->GetRenderWindow()->GetRenderers();
	
/*	std::cout << "rendercollection->GetNumberOfItems() = " << rendercollection->GetNumberOfItems() << std::endl;
	rendercollection->InitTraversal();
	for (vtkIdType i = 0; i < rendercollection->GetNumberOfItems(); i ++)
	{
	std::cout << "i = " << i << std::endl;
	vtkSmartPointer<vtkRenderer> thisrender = vtkRenderer::SafeDownCast(rendercollection->GetNextItem());
	vtkSmartPointer<vtkActorCollection> actorcollection = thisrender->GetActors();

	actorcollection->InitTraversal();
	for (vtkIdType j = 0; j < actorcollection->GetNumberOfItems(); j ++)
	{
	vtkSmartPointer<vtkActor> thisactor = vtkActor::SafeDownCast(actorcollection->GetNextActor());
	double color[3];
	thisactor->GetProperty()->GetColor(color);
	std::cout << i << ", " << j << ", " << color[0] << ", " << color[1] << ", " << color[2] << std::endl;
	}
	}
*/


	for (int i = 0; i < addedctrlobservertag.size(); i++)
		RenderWindowInteractorthreeD->RemoveObserver(addedctrlobservertag.at(i));
	addedctrlobservertag.clear();

	VesselPicker = vtkSmartPointer<vtkCellPicker>::New();
	VesselPicker->SetTolerance(0.005);
	VesselPicker->PickClippingPlanesOff();

	VesselPickCallBack = vtkSmartPointer<CVesselPickCallBack>::New();
	VesselPickCallBack->Slicer3DRenderWindowInteractor = RenderWindowInteractorthreeD;
	VesselPickCallBack->Slicer3DWidget = threeDView;
	VesselPickCallBack->clmodel = centerlineModel;
	VesselPickCallBack->lumenmodel = LumenModel;
	VesselPickCallBack->Slicer3DRender = rendercollection->GetFirstRenderer();
	VesselPickCallBack->logic = this;

	vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback> CtrlKeyPressedInteractionCallback = vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback>::New();
	CtrlKeyPressedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyPressedInteractionCallback->VesselPicker = VesselPicker;
	CtrlKeyPressedInteractionCallback->VesselPickCallBack = VesselPickCallBack;
	CtrlKeyPressedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	CtrlKeyPressedInteractionCallback->clmodel = centerlineModel;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyPressEvent, CtrlKeyPressedInteractionCallback));

	vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback> CtrlKeyReleasedInteractionCallback = vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback>::New();
	CtrlKeyReleasedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyReleasedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyReleaseEvent, CtrlKeyReleasedInteractionCallback));


	return true;
}



bool vtkSlicerCoronaryMainLogic
::TestLogic()
{
	std::cout << "TestLogic begin" << std::endl;

//	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
//	qMRMLThreeDView* threeDView = layoutManager->threeDWidget(0)->threeDView();
//	threeDView->VTKWidget()->setFocus();

	

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



















