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

#ifndef __qSlicerCoronaryMainModuleWidget_h
#define __qSlicerCoronaryMainModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerCoronaryMainModuleExport.h"

#include "qSlicerApplication.h"
#include "qSlicerLayoutManager.h"

// MRML includes
#include "vtkMRMLLinearTransformNode.h"
#include "vtkMRMLModelNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLNode.h"

#include "qMRMLLayoutManager.h"
#include "qMRMLThreeDWidget.h"
#include "qMRMLThreeDView.h"

// STD includes
#include <cstdlib>

// qt
#include "qsettings.h"
#include "qdir.h"
#include "qfiledialog.h"
#include "QShowEvent"

// VTK includes
#include <vtkObject.h>

#include "vtkXMLPolyDataReader.h"
#include "vtkMetaImageReader.h"

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

#include "vtkCellPicker.h"
#include "vtkSmartPointer.h"
#include "vtkMRMLScene.h"





class qSlicerCoronaryMainModuleWidgetPrivate;
class vtkMRMLNode;

class CVesselPickCallBack;

class QVesselEditingWidget : public QVTKWidget
{
	Q_OBJECT
public:
	QVesselEditingWidget(QVTKWidget *parent = 0, const char *name = 0){};
	~QVesselEditingWidget() {}

protected:
	virtual void mousePressEvent(QMouseEvent*);

public slots:
	void setvisibleslot(bool);
	
public:

};


/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_CORONARYMAIN_EXPORT qSlicerCoronaryMainModuleWidget :
  public qSlicerAbstractModuleWidget
{
	Q_OBJECT

public:

	typedef qSlicerAbstractModuleWidget Superclass;
	qSlicerCoronaryMainModuleWidget(QWidget *parent=0);
	virtual ~qSlicerCoronaryMainModuleWidget();

public:
	QString baseName;
	
public slots:
	bool DetectLandmarksButtonFunc();
	bool DetectCenterlinesButtonFunc();
	bool DetectLumenButtonFunc();
	bool SaveLandmarksButtonFunc();
	bool SaveCenterlinesButtonFunc();
	void SetVolumn(vtkMRMLNode* node);

	void SetCheckBoxBuildBifurcationMesh(int);

	void updateprogressbar(int i);

	bool TestButtonFunc();


public:
	void SavePolyData(vtkPolyData *poly, const char* fileName);
	void SaveVTKImage(vtkImageData *image, const char* fileName);

public:
	QVesselEditingWidget *VesselEditingWidget;

signals:
	void visibilitychanged(bool);

public:
	std::vector<vtkMRMLNode*> addedselectedclnode;
	std::vector< unsigned long > addedctrlobservertag;
	std::vector< unsigned long > addedvesselpickobservertag;
	vtkSmartPointer< vtkMRMLModelNode > SelectedClNode;
	vtkSmartPointer< vtkMRMLModelDisplayNode > SelectedClDisplayNode;
	vtkSmartPointer<vtkCellPicker> VesselPicker;
	vtkSmartPointer<CVesselPickCallBack> VesselPickCallBack;

public:
	bool RemoveAllSelectedVesselThreeD();
	bool ShowSelectedVesselThreeD(vtkIdType);
	bool SetupKeyMouseObserver();



protected:
	QScopedPointer<qSlicerCoronaryMainModuleWidgetPrivate> d_ptr;
	virtual void setup();

	vtkMRMLScalarVolumeNode* VolumeNode;
	vtkMRMLLinearTransformNode* TransformCoronaryNode;

private:
	Q_DECLARE_PRIVATE(qSlicerCoronaryMainModuleWidget);
	Q_DISABLE_COPY(qSlicerCoronaryMainModuleWidget);
};

#endif
