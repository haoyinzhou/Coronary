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

// MRML includes
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLModelNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLModelDisplayNode.h"

// VTK includes
#include <vtkObject.h>

// STD includes
#include <cstdlib>

#include "vtkSlicerCoronaryMainModuleLogicExport.h"

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


#include "Common.h"
#include "LearningImpl.h"
#include "ExtendTubeFilter.h"



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
  bool DetectLumenLogic();
  bool BuildMeshLogic();
  

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

private:

	Learning learn;
	vtkSmartPointer<vtkImageData> imageData;
	vtkSmartPointer<vtkImageData> imageData_original;
	vtkSmartPointer<vtkImageInterpolator> interpolator;
	double landmarks[SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS][3];
	vtkSmartPointer<vtkPolyData> centerlineModel;
	vtkSmartPointer<vtkPolyData> LumenModel;

	double NodeOrigin[3];
	double NodeSpaceing[3];


private:

  vtkSlicerCoronaryMainLogic(const vtkSlicerCoronaryMainLogic&); // Not implemented
  void operator=(const vtkSlicerCoronaryMainLogic&); // Not implemented
};

#endif
