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
#include "QDesktopWidget"
#include "QSpacerItem"


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
#include "QSizePolicy"

// VTK includes
#include <vtkObject.h>

#include "vtkSmartPointer.h"
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
#include "vtkLineSource.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageSlice.h"
#include "vtkImageProperty.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPropPicker.h"
#include "vtkAssemblyPath.h"
#include "vtkPointLocator.h"

#include "vtkCellPicker.h"
#include "vtkMRMLScene.h"
#include "ImageStretchCurvedReformat.h"
#include "ImageCurvedReformat.h"
#include "ImageObliqueReformat.h"
#include "Common.h"



class qSlicerCoronaryMainModuleWidgetPrivate;
class vtkMRMLNode;

class CVesselPickCallBack;
class ORSliceStyle;

class QVesselEditingWidget : public QVTKWidget
{
	Q_OBJECT
public:
	QVesselEditingWidget();
	~QVesselEditingWidget();

public:
	QVTKWidget* widget1;
	QVTKWidget* widget2;

	vtkSmartPointer<ORSliceStyle> ORSliceStyleCallback;


signals:
	void clcoordchanged(vtkIdType, double, double, double);
	void lumenradiuschanged(vtkIdType, vtkIdType, double);
public:
	void send_clcoordchanged(vtkIdType, double, double, double);
	void send_lumenradiuschanged(vtkIdType, vtkIdType, double);

	
public slots:
	void setvisibleslot(bool);
	void setselectidslot(vtkIdType);
	void setclmodelslot(vtkPolyData*);
	void setimagedataslot(vtkImageData*);
	void resetslot(void);
	void forcerenderslot(void);

public:
	vtkPolyData* clModel;
	vtkImageData* ImageData;
	vtkIdType SelectID;

	vtkSmartPointer<ImageStretchCurvedReformat>	 stretchCurvedReformat;
	vtkSmartPointer<ImageCurvedReformat>	 CurvedReformat;
	vtkSmartPointer<ImageObliqueReformat> ObliqueReformat;

public:
	void SavePolyData(vtkPolyData *poly, const char* fileName);
	void SaveVTKImage(vtkImageData *image, const char* fileName);  // just for debug
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
	void selectidchanged(vtkIdType);
	void clmodelchanged(vtkPolyData*);
	void imagedatachanged(vtkImageData*);
	void resetsingal();
	void forcerendersingal();

public:
	void send_visibilitychanged(bool);
	void send_selectidchanged(vtkIdType);
	void send_clmodelchanged(vtkPolyData*);
	void send_imagedatachanged(vtkImageData*);
	void send_resetsingal();
	void send_forcerendersingal();

public slots:
	void setclcoordslot(vtkIdType, double, double, double);
	void setlumenradiusslot(vtkIdType, vtkIdType, double);


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
