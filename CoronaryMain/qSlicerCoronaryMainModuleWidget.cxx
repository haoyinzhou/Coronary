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

// Qt includes
#include <QDebug>
#include <QMessageBox>

// SlicerQt includes
#include "qSlicerCoronaryMainModuleWidget.h"
#include "ui_qSlicerCoronaryMainModuleWidget.h"
#include "vtkSlicerCoronaryMainLogic.h"




//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerCoronaryMainModuleWidgetPrivate: public Ui_qSlicerCoronaryMainModuleWidget
{
	Q_DECLARE_PUBLIC(qSlicerCoronaryMainModuleWidget);

protected:
	qSlicerCoronaryMainModuleWidget* const q_ptr;

public:
	qSlicerCoronaryMainModuleWidgetPrivate(qSlicerCoronaryMainModuleWidget& object);
	vtkSlicerCoronaryMainLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerCoronaryMainModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModuleWidgetPrivate::qSlicerCoronaryMainModuleWidgetPrivate(qSlicerCoronaryMainModuleWidget& object) : q_ptr(&object)
{
}

vtkSlicerCoronaryMainLogic* qSlicerCoronaryMainModuleWidgetPrivate::logic() const
{
	Q_Q(const qSlicerCoronaryMainModuleWidget);
	return vtkSlicerCoronaryMainLogic::SafeDownCast(q->logic());
}
//-----------------------------------------------------------------------------
// qSlicerCoronaryMainModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModuleWidget::qSlicerCoronaryMainModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerCoronaryMainModuleWidgetPrivate(*this) )
{
	
}

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModuleWidget::~qSlicerCoronaryMainModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerCoronaryMainModuleWidget::setup()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	d->setupUi(this);
	this->Superclass::setup();

	this->VolumeNode = NULL;
	this->TransformCoronaryNode = NULL;

	connect(d->DetectLandmarks, SIGNAL(clicked()), this, SLOT(DetectLandmarksButtonFunc()));
	connect(d->DetectCenterlines, SIGNAL(clicked()), this, SLOT(DetectCenterlinesButtonFunc()));
	connect(d->DetectLumen, SIGNAL(clicked()), this, SLOT(DetectLumenButtonFunc()));
	connect(d->SaveLandmarks, SIGNAL(clicked()), this, SLOT(SaveLandmarksButtonFunc()));
	connect(d->SaveCenterlines, SIGNAL(clicked()), this, SLOT(SaveCenterlinesButtonFunc()));
	connect(d->MRMLNodeReadVolumn, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(SetVolumn(vtkMRMLNode*)));
	connect(d->checkBox_buildbifurcationmesh, SIGNAL(stateChanged(int)), this, SLOT(SetCheckBoxBuildBifurcationMesh(int)));

	d->progressBar->setValue(0);
}

#define DisableAllButtons()\
{\
	d->DetectLandmarks->setEnabled(false);\
	d->DetectCenterlines->setEnabled(false);\
	d->DetectLumen->setEnabled(false);\
}\

#define EnableAllButtons()\
{\
	d->DetectLandmarks->setEnabled(true);\
	d->DetectCenterlines->setEnabled(true);\
	d->DetectLumen->setEnabled(true);\
}\


void qSlicerCoronaryMainModuleWidget::SetVolumn(vtkMRMLNode* node)
{
	if (node != NULL)
	{
		this->VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(node);
	}
}

void qSlicerCoronaryMainModuleWidget::SetCheckBoxBuildBifurcationMesh(int state)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic == NULL)
		return;

	if (state)
	{
		std::cout << "Will build bifurcation mesh" << std::endl;
		logic->WillBuildBifurcationMesh = true;
	}
	else
	{
		std::cout << "Will not build bifurcation mesh" << std::endl;
		logic->WillBuildBifurcationMesh = false;
	}

	logic->BuildMeshLogic();
}


bool qSlicerCoronaryMainModuleWidget::DetectLandmarksButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
	//	DisableAllButtons();
		logic->DetectLandmarksLogic(VolumeNode, d->progressBar);
	//	EnableAllButtons();
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::DetectCenterlinesButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		if (logic->DetectCenterlinesLogic(d->progressBar))
			logic->BuildMeshLogic();
	}
	return true;
}

bool qSlicerCoronaryMainModuleWidget::DetectLumenButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		logic->DetectLumenLogic();
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::SaveLandmarksButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		logic->SaveLandmarksLogic();
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::SaveCenterlinesButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		logic->SaveCenterlinesLogic();
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::BuildMesh()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		logic->BuildMeshLogic();
	}
	return true;
}

void qSlicerCoronaryMainModuleWidget::updateprogressbar(int i)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	d->progressBar->setValue(i);
}

