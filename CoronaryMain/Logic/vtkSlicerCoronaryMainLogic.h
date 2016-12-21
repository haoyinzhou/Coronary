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

// .NAME vtkSlicerCoronaryMainLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerCoronaryMainLogic_h
#define __vtkSlicerCoronaryMainLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// SlicerQt includes
#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"

// MRML includes
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLModelNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLModelDisplayNode.h"

#include "vtkMRMLCrosshairNode.h"

#include "qMRMLLayoutManager.h"
#include "qMRMLThreeDWidget.h"
#include "qMRMLThreeDView.h"


// VTK includes
#include <vtkObject.h>

// STD includes
#include <cstdlib>

#include "vtkSlicerCoronaryMainModuleLogicExport.h"
#include "qSlicerAbstractModuleWidget.h"

#include "vtkImageInterpolator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkMetaImageWriter.h"
#include "vtkSmartPointer.h"
#include "vtkSphereWidget.h"
#include "vtkCaptionActor2D.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkCollection.h"

#include "vtkCornerAnnotation.h"
#include "vtkTextProperty.h"

#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkPointPicker.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkObjectFactory.h"

#include "vtkSelectionNode.h"
#include "vtkSelection.h"
#include "vtkExtractSelection.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGeometryFilter.h"
#include "vtkCellPicker.h"
#include "vtkRendererCollection.h"
#include "vtkProperty.h"

#include "Common.h"
#include "LearningImpl.h"
#include "ExtendTubeFilter.h"


class CVesselPickCallBack;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_CORONARYMAIN_MODULE_LOGIC_EXPORT vtkSlicerCoronaryMainLogic :
  public vtkSlicerModuleLogic
{
public:
  static vtkSlicerCoronaryMainLogic *New();
  vtkTypeMacro(vtkSlicerCoronaryMainLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  bool DetectLandmarksLogic(vtkMRMLScalarVolumeNode* VolumnNode, QProgressBar* progressbar);
  bool DetectCenterlinesLogic(QProgressBar* progressbar);
  bool DetectLumenLogic(QProgressBar* progressbar);

  bool BuildLandmarksMeshLogic();
  bool BuildCenterlinesMeshLogic();

  bool SetupKeyMouseObserver();
  bool ShowSelectedVesselThreeD(vtkIdType);
  bool RemoveAllSelectedVesselThreeD();
  
  bool TestLogic();  // just for debug

public:
	bool GetLandMarksCoord(int index, double coord[3]);
	bool SetLandMarksCoord(int index, double coord[3]);

protected:
  vtkSlicerCoronaryMainLogic();
  virtual ~vtkSlicerCoronaryMainLogic();

  /// Initialize listening to MRML events
  virtual void SetMRMLSceneInternal(vtkMRMLScene * newScene);
  virtual void ObserveMRMLScene();

  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);


public:	
	Learning learn;
	vtkSmartPointer<vtkImageData> imageData;
	vtkSmartPointer<vtkImageData> imageData_original;
	vtkSmartPointer<vtkImageData> hessianImage;
	vtkSmartPointer<vtkImageInterpolator> interpolator;
	double landmarks[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS][3];
	vtkSmartPointer<vtkPolyData> centerlineModel;
	vtkSmartPointer<vtkPolyData> LumenModel;
	vtkSmartPointer<vtkPolyData> centerlineModel_display;
	vtkSmartPointer<vtkPolyData> LumenModel_display;

	vtkSmartPointer<vtkIdFilter> centerlineId;
	vtkSmartPointer<ExtendTubeFilter> centerlineTube;

	double NodeOrigin[3];
	double NodeSpaceing[3];
	
	vector<vtkMRMLNode*> addedclnode;
	vector<vtkMRMLNode*> addedlandmarknode;
	vector<vtkMRMLNode*> addedselectedclnode;

	vector< unsigned long > addedctrlobservertag;
	vector< unsigned long > addedvesselpickobservertag;

	vtkSmartPointer< vtkMRMLModelNode > LandmarkNode;
	vtkSmartPointer< vtkMRMLModelDisplayNode > LandmarkDisplayNode;
	vtkSmartPointer< vtkMRMLModelNode > clNode;
	vtkSmartPointer< vtkMRMLModelDisplayNode > clDisplayNode;
	vtkSmartPointer< vtkMRMLModelNode > LumenNode;
	vtkSmartPointer< vtkMRMLModelDisplayNode > LumenDisplayNode;

	vtkSmartPointer<vtkRenderer> LandmarksRender;
	vtkSmartPointer<vtkRenderer> clRender;
	vtkSmartPointer<vtkRenderer> LumenRender;

	vtkSmartPointer< vtkMRMLModelNode > SelectedClNode;
	vtkSmartPointer< vtkMRMLModelDisplayNode > SelectedClDisplayNode;

	vtkSmartPointer<vtkCellPicker> VesselPicker;
	vtkSmartPointer<CVesselPickCallBack> VesselPickCallBack;

private:

public:
	bool WillBuildBifurcationMesh;


public:
	vtkIdType cellid_temp; // just for debug
	
private:

  vtkSlicerCoronaryMainLogic(const vtkSlicerCoronaryMainLogic&); // Not implemented
  void operator=(const vtkSlicerCoronaryMainLogic&); // Not implemented
};






#endif
