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



QVesselEditingWidget::QVesselEditingWidget()
{

}

QVesselEditingWidget::~QVesselEditingWidget()
{

}


void QVesselEditingWidget::mousePressEvent(QMouseEvent *	e)
{
	if (e->button() == Qt::LeftButton)
	{

	}
	QVTKWidget::mousePressEvent(e);
}

void QVesselEditingWidget::setvisibleslot(bool f)
{
	std::cout << "set visible slot" << std::endl;
	this->setVisible(f);
}

void QVesselEditingWidget::setselectidslot(vtkIdType id)
{
	std::cout << "set selectid slot, id = " << id << std::endl;
	this->SelectID = id;
}

void QVesselEditingWidget::setclmodelslot(vtkPolyData* cl, vtkPolyData* lumen)
{
	std::cout << "set clmodel slot" << std::endl;
	this->clModel = cl;
	this->lumenModel = lumen;
}

void QVesselEditingWidget::setimagedataslot(vtkImageData* p)
{
	std::cout << "set imagedata slot" << std::endl;
	this->ImageData = p;
}

void QVesselEditingWidget::resetslot()
{
	std::cout << "reset slot" << std::endl;

	this->clModel = NULL;
	this->ImageData = NULL;

	vtkSmartPointer<vtkRendererCollection> rendercollection = this->GetInteractor()->GetRenderWindow()->GetRenderers();
	
	rendercollection->InitTraversal();
	for (vtkIdType i = 0; i < rendercollection->GetNumberOfItems(); i++)
	{
		vtkSmartPointer<vtkRenderer> thisrender = vtkRenderer::SafeDownCast(rendercollection->GetNextItem());
		this->GetInteractor()->GetRenderWindow()->RemoveRenderer(thisrender);
	}

}

void QVesselEditingWidget::forcerenderslot()
{
	std::cout << "force render slot" << std::endl;
//	SaveVTKImage(this->ImageData, "C:\\work\\Coronary_Slicer\\testdata\\imageinforcerenderslot.mha");


	stretchCurvedReformat = vtkSmartPointer<ImageStretchCurvedReformat>::New();
	stretchCurvedReformatLine = vtkSmartPointer<vtkLineSource>::New();
	vtkSmartPointer<vtkActor> stretchCurvedReformatActor = vtkSmartPointer<vtkActor>::New();

//	SaveVTKImage(ImageData, "C:\\work\\Coronary_Slicer\\testdata\\image_input.mha");
//	SavePolyData(clModel, "C:\\work\\Coronary_Slicer\\testdata\\clModel_input.vtp");

	std::cout << "SelectID = " << SelectID << std::endl;


	stretchCurvedReformat->SetInputData(0, ImageData);
	stretchCurvedReformat->SetInputData(1, clModel);
	stretchCurvedReformat->SetSegmentId(&SelectID);
	stretchCurvedReformat->SetTwistIndex(0);
	stretchCurvedReformat->Update();

	
	SaveVTKImage(stretchCurvedReformat->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output0.mha");
	vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(stretchCurvedReformat->GetOutput(1));
	SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output1.vtp");
	poly = vtkPolyData::SafeDownCast(stretchCurvedReformat->GetOutput(2));
	SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output2.vtp");
	poly = vtkPolyData::SafeDownCast(stretchCurvedReformat->GetOutput(3));
	SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output3.vtp");
	poly = vtkPolyData::SafeDownCast(stretchCurvedReformat->GetOutput(4));
	SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output4.vtp");
	poly = vtkPolyData::SafeDownCast(stretchCurvedReformat->GetOutput(5));
	SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\stretchCurvedReformat_output5.vtp");



/*	vtkSmartPointer<vtkPolyDataMapper> circleContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	circleContourMapper->SetInputConnection(stretchCurvedReformat->GetOutputPort(1));
	vtkSmartPointer<vtkActor> circleContourActor = vtkSmartPointer<vtkActor>::New();
	circleContourActor->SetMapper(circleContourMapper);
	circleContourActor->GetProperty()->SetColor(0.8, 0.0, 0.0);
	circleContourActor->GetProperty()->SetLineWidth(2.0f);
	circleContourActor->GetProperty()->SetOpacity(0.6);
	circleContourActor->PickableOff();
	circleContourActor->VisibilityOff();
	vtkSmartPointer<vtkRenderer> render1 = vtkSmartPointer<vtkRenderer>::New();
	render1->AddActor(circleContourActor);
	this->GetRenderWindow()->AddRenderer(render1);
*/




/*	vtkSmartPointer<vtkPolyDataMapper> VesselEditingMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	VesselEditingMapper->SetInputData(lumenModel);
	vtkSmartPointer<vtkActor> VesselEditingActor = vtkSmartPointer<vtkActor>::New();
	VesselEditingActor->SetMapper(VesselEditingMapper);
	vtkSmartPointer<vtkRenderer> VesselEditingRenderer = vtkSmartPointer<vtkRenderer>::New();
	VesselEditingRenderer->SetBackground(0.0, 0.0, 0.0); // Background color black
	VesselEditingRenderer->AddActor(VesselEditingActor);

	this->GetRenderWindow()->AddRenderer(VesselEditingRenderer);
*/


	//vtkSmartPointer<vtkRenderWindowInteractor> VesselEditingRenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//VesselEditingRenderWindowInteractor->SetRenderWindow(VesselEditingRenderWindow);

	//vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	//renderWindowInteractor->SetInteractorStyle(style);

	//VesselEditingRenderWindowInteractor->Start();

}

void QVesselEditingWidget::SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}


void QVesselEditingWidget::SaveVTKImage(vtkImageData *image, const char* fileName)
{
	vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
	writer->SetFileName(fileName);
	writer->SetInputData(image);
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}





void qSlicerCoronaryMainModuleWidget::send_visibilitychanged(bool f)
{
	emit visibilitychanged(f);
}
void qSlicerCoronaryMainModuleWidget::send_selectidchanged(vtkIdType id)
{
	emit selectidchanged(id);
}
void qSlicerCoronaryMainModuleWidget::send_clmodelchanged(vtkPolyData* cl, vtkPolyData* lumen)
{
	emit clmodelchanged(cl, lumen);
}
void qSlicerCoronaryMainModuleWidget::send_imagedatachanged(vtkImageData* p)
{
	emit imagedatachanged(p);
}
void qSlicerCoronaryMainModuleWidget::send_resetsingal()
{
	emit resetsingal();
}
void qSlicerCoronaryMainModuleWidget::send_forcerendersingal()
{
	emit forcerendersingal();
}


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
	baseName = "dataset";
}

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModuleWidget::~qSlicerCoronaryMainModuleWidget()
{
	VesselEditingWidget->deleteLater();
	delete[] VesselEditingWidget;
	baseName.clear();
}

//-----------------------------------------------------------------------------
void qSlicerCoronaryMainModuleWidget::setup()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	d->setupUi(this);
	this->Superclass::setup();

	this->VolumeNode = NULL;
	this->TransformCoronaryNode = NULL;

	VesselEditingWidget = new QVesselEditingWidget;
	VesselEditingWidget->resize(600, 1200);
	VesselEditingWidget->setVisible(false);


	connect(d->DetectLandmarks, SIGNAL(clicked()), this, SLOT(DetectLandmarksButtonFunc()));
	connect(d->DetectCenterlines, SIGNAL(clicked()), this, SLOT(DetectCenterlinesButtonFunc()));
	connect(d->DetectLumen, SIGNAL(clicked()), this, SLOT(DetectLumenButtonFunc()));
	connect(d->SaveLandmarks, SIGNAL(clicked()), this, SLOT(SaveLandmarksButtonFunc()));
	connect(d->SaveCenterlines, SIGNAL(clicked()), this, SLOT(SaveCenterlinesButtonFunc()));
	connect(d->MRMLNodeReadVolumn, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(SetVolumn(vtkMRMLNode*)));
	connect(d->checkBox_buildbifurcationmesh, SIGNAL(stateChanged(int)), this, SLOT(SetCheckBoxBuildBifurcationMesh(int)));

	connect(this, SIGNAL(visibilitychanged(bool)), VesselEditingWidget, SLOT(setvisibleslot(bool)));
	connect(this, SIGNAL(selectidchanged(vtkIdType)), VesselEditingWidget, SLOT(setselectidslot(vtkIdType)));
	connect(this, SIGNAL(clmodelchanged(vtkPolyData*, vtkPolyData*)), VesselEditingWidget, SLOT(setclmodelslot(vtkPolyData*, vtkPolyData*)));
	connect(this, SIGNAL(imagedatachanged(vtkImageData*)), VesselEditingWidget, SLOT(setimagedataslot(vtkImageData*)));
	connect(this, SIGNAL(resetsingal(void)), VesselEditingWidget, SLOT(resetslot()));
	connect(this, SIGNAL(forcerendersingal(void)), VesselEditingWidget, SLOT(forcerenderslot()));

	
	connect(d->pushButtonTest, SIGNAL(clicked()), this, SLOT(TestButtonFunc()));

	d->progressBar->setValue(0);
	d->checkBox_buildbifurcationmesh->setChecked(false);
	d->checkBox_loadlandmarks->setChecked(false);
	d->checkBox_loadcenterlines->setChecked(false);

	addedselectedclnode.clear();
	addedctrlobservertag.clear();
	addedvesselpickobservertag.clear();
}


void qSlicerCoronaryMainModuleWidget::SetVolumn(vtkMRMLNode* node)
{
	if (node != NULL)
	{
		Q_D(qSlicerCoronaryMainModuleWidget);
		vtkSlicerCoronaryMainLogic *logic = d->logic();

		this->VolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(node);

		this->VolumeNode->GetOrigin(logic->NodeOrigin);
		this->VolumeNode->GetSpacing(logic->NodeSpaceing);

		logic->imageData_original = this->VolumeNode->GetImageData();
		logic->imageData->DeepCopy(logic->imageData_original);
		logic->imageData->SetOrigin(-(logic->NodeOrigin[0]), -(logic->NodeOrigin[1]), logic->NodeOrigin[2]);
		logic->imageData->SetSpacing(logic->NodeSpaceing);

		logic->interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
		logic->interpolator->SetInterpolationModeToLinear();
		logic->interpolator->SetOutValue(-3024.0);
		logic->interpolator->Initialize(logic->imageData);
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
		logic->WillBuildBifurcationMesh = true;
	}
	else
	{
		logic->WillBuildBifurcationMesh = false;
	}

	if (logic->centerlineModel->GetNumberOfCells() == 0)
		return;

	logic->BuildCenterlinesMeshLogic();
	SetupKeyMouseObserver();

}


bool qSlicerCoronaryMainModuleWidget::DetectLandmarksButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		if (d->checkBox_loadlandmarks->isChecked() == false)
		{
			logic->DetectLandmarksLogic(VolumeNode, d->progressBar);
			for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
			{
				std::cout << "landmarkcoord " << i << " = " << logic->landmarks[i][0] << ", " << logic->landmarks[i][1] << ", " << logic->landmarks[i][2] << std::endl;
			}
		}
		else
		{
			std::cout << "load landmarks..." << std::endl;
			d->progressBar->setValue(0);

			const QString DEFAULT_DIR_KEY("default_dir");
			QSettings MySettings;
			QString suggestName(QDir::separator());
		//	suggestName += baseName + tr("-lvcorlm.vtp");
			suggestName += tr("landmarks.vtp");
			QString fileName = QFileDialog::getOpenFileName(this, tr("Open Landmark File"), MySettings.value(DEFAULT_DIR_KEY).toString() + suggestName, tr("Landmark File (*.vtp)"));
			if (fileName.isEmpty()) return false;
			QByteArray fileNameByte = fileName.toLocal8Bit();
			vtkSmartPointer< vtkXMLPolyDataReader > reader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
			reader->SetFileName(fileNameByte.data());
			try
			{
				reader->Update();
			}
			catch (...)
			{
				std::cerr << "Error occurs when reading " << fileNameByte.data() << std::endl;
				return false;
			}

			QFileInfo fileInfo(fileName);
		//	baseName = fileInfo.baseName();
			MySettings.setValue(DEFAULT_DIR_KEY, fileInfo.absolutePath());

			if (reader->GetOutput()->GetNumberOfPoints() < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS)
			{
				std::cerr << "The number of points in the file " << fileNameByte.data() << " is incorrect" << std::endl;
				return false;
			}
			double coord[3];
			for (int i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
			{
				reader->GetOutput()->GetPoint(i, coord);
				std::cout << "landmarkcoord " << i << " = " << coord[0] << ", " << coord[1] << ", " << coord[2] << std::endl;
				logic->SetLandMarksCoord(i, coord);
			}
			d->progressBar->setValue(100);
		}
		logic->BuildLandmarksMeshLogic();
	}
	
	return true;
}

bool qSlicerCoronaryMainModuleWidget::DetectCenterlinesButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (logic != NULL)
	{
		if (logic->imageData->GetNumberOfCells() == 0)
		{
			std::cerr << "cannot find image data" << std::endl;
			return false;
		}

		if (d->checkBox_loadcenterlines->isChecked() == false)
		{
			double z[3] = { 0.0, 0.0, 0.0 };
			if (vtkMath::Distance2BetweenPoints(logic->landmarks[SmartCoronary::LEFT_CORONARY_OSTIUM], z) < 1e-3
				&& vtkMath::Distance2BetweenPoints(logic->landmarks[SmartCoronary::RIGHT_CORONARY_OSTIUM], z) < 1e-3)
			{
				std::cerr << "cannot find landmark!" << std::endl;
				return false;
			}
			logic->centerlineModel = vtkSmartPointer<vtkPolyData>::New();
			logic->centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel_display = vtkSmartPointer<vtkPolyData>::New();
			logic->DetectCenterlinesLogic(d->progressBar);
		}
		else
		{
			logic->centerlineModel = vtkSmartPointer<vtkPolyData>::New();
			logic->centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel_display = vtkSmartPointer<vtkPolyData>::New();

			d->progressBar->setValue(0);
			const QString DEFAULT_DIR_KEY("default_dir");
			QSettings MySettings;
			QString suggestName(QDir::separator());
		//	suggestName += baseName + tr("-centerline.vtp");
			suggestName += tr("centerlines.vtp");
			QString fileName = QFileDialog::getOpenFileName(this, tr("Open Centerline File"), MySettings.value(DEFAULT_DIR_KEY).toString() + suggestName, tr("Landmark File (*.vtp)"));
			if (fileName.isEmpty()) return false;
			QByteArray fileNameByte = fileName.toLocal8Bit();
			vtkSmartPointer< vtkXMLPolyDataReader > reader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
			reader->SetFileName(fileNameByte.data());
			try
			{
				reader->Update();
			}
			catch (...)
			{
				std::cerr << "Error occurs when reading " << fileNameByte.data() << std::endl;
				return false;
			}

			QFileInfo fileInfo(fileName);
		//	baseName = fileInfo.baseName();
			MySettings.setValue(DEFAULT_DIR_KEY, fileInfo.absolutePath());
		//	std::cout << "baseName = " << baseName.toStdString() << endl;

			logic->centerlineModel = vtkSmartPointer<vtkPolyData>::New();
			logic->centerlineModel->DeepCopy(reader->GetOutput());
			//remove "SegmentId" if there is one in the loaded centerline because we will generate the ids with idfilter
			if (logic->centerlineModel->GetCellData()->HasArray("SegmentId"))	logic->centerlineModel->GetCellData()->RemoveArray("SegmentId");

			QFileInfo fi(fileName);
			//	QString fileName2 = fi.absolutePath() + "/" + fi.baseName() + "-hessian.mha";
			QString fileName2 = fi.absolutePath() + "/" + "hessian.mha";
			QFile hessianFile(fileName2);
			if (hessianFile.exists())
			{
				fileNameByte = fileName2.toLocal8Bit();
				vtkSmartPointer< vtkMetaImageReader > hessianReader = vtkSmartPointer< vtkMetaImageReader >::New();
				hessianReader->SetFileName(fileNameByte.data());
				try
				{
					hessianReader->Update();
					logic->hessianImage = vtkSmartPointer<vtkImageData>::New();
					logic->hessianImage->DeepCopy(hessianReader->GetOutput());
				}
				catch (...)
				{
					std::cerr << "Error occurs when reading " << fileNameByte.data() << std::endl;
				}
			}
			else
			{
				std::cerr << "cannot find hessianImage" << std::endl;
				return false;
			}

			d->progressBar->setValue(95);
		}

		logic->centerlineId = vtkSmartPointer<vtkIdFilter>::New();
		logic->centerlineId->SetInputData(logic->centerlineModel);
		logic->centerlineId->PointIdsOff();
		logic->centerlineId->CellIdsOn();
		logic->centerlineId->FieldDataOn();
		logic->centerlineId->SetIdsArrayName("SegmentId");
		logic->centerlineId->Update();
		logic->centerlineModel = vtkPolyData::SafeDownCast(logic->centerlineId->GetOutput());
		
		logic->BuildCenterlinesMeshLogic();
				
		SetupKeyMouseObserver();

		d->progressBar->setValue(100);
	}
	return true;
}

bool qSlicerCoronaryMainModuleWidget::DetectLumenButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		if (logic->DetectLumenLogic(d->progressBar))
		{
			logic->BuildCenterlinesMeshLogic();
			SetupKeyMouseObserver();
		}
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::SaveLandmarksButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic != NULL)
	{
		const QString DEFAULT_DIR_KEY("default_dir");
		QSettings MySettings;
		QString suggestName(QDir::separator());
	//	suggestName += baseName + tr("-lvcorlm.vtp");
		suggestName += tr("landmarks.vtp");
		QString fileName = QFileDialog::getSaveFileName(this, tr("Save Landmark File"), MySettings.value(DEFAULT_DIR_KEY).toString() + suggestName, tr("Model Files (*.vtp)"));
		if (fileName.isEmpty()) return false;
		QByteArray fileNameByte = fileName.toLocal8Bit();

		QFileInfo fileInfo(fileName);
	//	baseName = fileInfo.baseName();
		MySettings.setValue(DEFAULT_DIR_KEY, fileInfo.absolutePath());

		vtkSmartPointer<vtkPolyData> lmPoly = vtkSmartPointer<vtkPolyData>::New();
		lmPoly->Allocate(SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS);
		vtkSmartPointer<vtkPoints> lmPoints = vtkSmartPointer<vtkPoints>::New();
		lmPoints->SetNumberOfPoints(SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS);
		double landmarkcoord[3];
		for (vtkIdType i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
		{
			logic->GetLandMarksCoord(i, landmarkcoord);
			std::cout << "landmarkcoord " << i << " = " << landmarkcoord[0] << ", " << landmarkcoord[1] << ", " << landmarkcoord[2] << std::endl;
			lmPoints->SetPoint(i, landmarkcoord);
		}
		for (vtkIdType i = 0; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
			lmPoly->InsertNextCell(VTK_VERTEX, 1, &i);
		lmPoly->SetPoints(lmPoints);
		std::cout << "save " << fileNameByte.data() << std::endl;
		qSlicerCoronaryMainModuleWidget::SavePolyData(lmPoly, fileNameByte.data());
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget::SaveCenterlinesButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	if (logic == NULL)
		return false;
	if (logic->centerlineModel->GetPointData()->GetNumberOfArrays() == 0)
	{
		std::cerr << "cannot find centerline model" << std::endl;
		return false;
	}	

	const QString DEFAULT_DIR_KEY("default_dir");
	QSettings MySettings;
	QString suggestName(QDir::separator());
	//suggestName += baseName + tr("-centerline.vtp");
	suggestName += tr("centerlines.vtp");
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Centerline File"), MySettings.value(DEFAULT_DIR_KEY).toString() + suggestName, tr("Model Files (*.vtp)"));
	if (fileName.isEmpty()) return false;
	QByteArray fileNameByte = fileName.toLocal8Bit();

	QFileInfo fileInfo(fileName);
//	baseName = fileInfo.baseName();
	MySettings.setValue(DEFAULT_DIR_KEY, fileInfo.absolutePath());

	SavePolyData(logic->centerlineModel, fileNameByte.data());

	if (logic->hessianImage)
	{
		QFileInfo fi(fileName);
		//	QString fileName2 = fi.absolutePath() + "/" + fi.baseName() + "-hessian.mha";
		QString fileName2 = fi.absolutePath() + "/" + "hessian.mha";
		QByteArray fileNameByte = fileName2.toLocal8Bit();
		qSlicerCoronaryMainModuleWidget::SaveVTKImage(logic->hessianImage, fileNameByte.data());
	}

	if (logic->LumenModel->GetPointData()->GetNumberOfArrays() != 0)
	{
		QFileInfo fi(fileName);

		QString fileName2;
		{
		//	fileName2 = fi.absolutePath() + "/" + fi.baseName() + "-LumenModel.vtp";
			fileName2 = fi.absolutePath() + "/" + "LumenModel.vtp";
			QByteArray fileNameByte = fileName2.toLocal8Bit();
			SavePolyData(logic->LumenModel, fileNameByte.data());
		}
	}

	return true;
}



void qSlicerCoronaryMainModuleWidget::updateprogressbar(int i)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	d->progressBar->setValue(i);
}

void qSlicerCoronaryMainModuleWidget::SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

void qSlicerCoronaryMainModuleWidget::SaveVTKImage(vtkImageData *image, const char* fileName)
{
	vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
	writer->SetFileName(fileName);
	writer->SetInputData(image);
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}


bool qSlicerCoronaryMainModuleWidget
::RemoveAllSelectedVesselThreeD()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (addedselectedclnode.size() != 0)
	{
		for (int i = 0; i < addedselectedclnode.size(); i++)
			logic->GetMRMLScene()->RemoveNode(addedselectedclnode.at(i));
		addedselectedclnode.clear();
	}

	return true;
}

bool qSlicerCoronaryMainModuleWidget
::ShowSelectedVesselThreeD(vtkIdType cellid)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	std::cout << "ShowSelectedVesselThreeD begin!" << std::endl;

	vtkSmartPointer<vtkIdTypeArray> centerlineSelectId = vtkSmartPointer<vtkIdTypeArray>::New();
	centerlineSelectId->SetName("SegmentId");
	centerlineSelectId->SetNumberOfValues(1);
	centerlineSelectId->SetValue(0, cellid);

	if (cellid >= logic->centerlineModel->GetNumberOfCells())
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
	extractSelection->SetInputData(0, logic->LumenModel_display);
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

	thisaddednode = logic->GetMRMLScene()->AddNode(SelectedClNode);
	addedselectedclnode.push_back(thisaddednode);
	thisaddednode = logic->GetMRMLScene()->AddNode(SelectedClDisplayNode);
	addedselectedclnode.push_back(thisaddednode);
	SelectedClDisplayNode->SetScene(logic->GetMRMLScene());
	SelectedClNode->SetScene(logic->GetMRMLScene());
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

		mainwidget->ShowSelectedVesselThreeD(pickid);


		mainwidget->send_visibilitychanged(true);
		mainwidget->send_selectidchanged(pickid);
		mainwidget->send_clmodelchanged(clmodel, lumenmodel);
		mainwidget->send_imagedatachanged(imagedata);
		//mainwidget->send_resetsingal();
		mainwidget->send_forcerendersingal();

	}

public:
	vtkRenderWindowInteractor* Slicer3DRenderWindowInteractor;
	vtkRenderer* Slicer3DRender;
	qSlicerCoronaryMainModuleWidget* mainwidget;
	
	vtkPolyData* clmodel;
	vtkPolyData* lumenmodel;
	vtkImageData* imagedata;

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



bool qSlicerCoronaryMainModuleWidget
::SetupKeyMouseObserver()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (logic->centerlineModel->GetNumberOfCells() == 0)
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
	VesselPickCallBack->Slicer3DRender = rendercollection->GetFirstRenderer();
	VesselPickCallBack->mainwidget = this;
	VesselPickCallBack->clmodel = logic->centerlineModel;
	VesselPickCallBack->lumenmodel = logic->LumenModel;
	VesselPickCallBack->imagedata = logic->imageData;


	vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback> CtrlKeyPressedInteractionCallback = vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback>::New();
	CtrlKeyPressedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyPressedInteractionCallback->VesselPicker = VesselPicker;
	CtrlKeyPressedInteractionCallback->VesselPickCallBack = VesselPickCallBack;
	CtrlKeyPressedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	CtrlKeyPressedInteractionCallback->clmodel = logic->centerlineModel;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyPressEvent, CtrlKeyPressedInteractionCallback));

	vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback> CtrlKeyReleasedInteractionCallback = vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback>::New();
	CtrlKeyReleasedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyReleasedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyReleaseEvent, CtrlKeyReleasedInteractionCallback));


	return true;
}




bool qSlicerCoronaryMainModuleWidget::TestButtonFunc()
{
	std::cout << "TestButtonFunc begin" << std::endl;

	Q_D(qSlicerCoronaryMainModuleWidget);

	d->logic()->TestLogic();

	std::cout << "TestButtonFunc end" << std::endl;

	return true;
}
