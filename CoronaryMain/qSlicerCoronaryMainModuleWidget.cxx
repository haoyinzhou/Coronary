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
	baseName = "dataset";
}

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModuleWidget::~qSlicerCoronaryMainModuleWidget()
{
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

	connect(d->DetectLandmarks, SIGNAL(clicked()), this, SLOT(DetectLandmarksButtonFunc()));
	connect(d->DetectCenterlines, SIGNAL(clicked()), this, SLOT(DetectCenterlinesButtonFunc()));
	connect(d->DetectLumen, SIGNAL(clicked()), this, SLOT(DetectLumenButtonFunc()));
	connect(d->SaveLandmarks, SIGNAL(clicked()), this, SLOT(SaveLandmarksButtonFunc()));
	connect(d->SaveCenterlines, SIGNAL(clicked()), this, SLOT(SaveCenterlinesButtonFunc()));
	connect(d->MRMLNodeReadVolumn, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(SetVolumn(vtkMRMLNode*)));
	connect(d->checkBox_buildbifurcationmesh, SIGNAL(stateChanged(int)), this, SLOT(SetCheckBoxBuildBifurcationMesh(int)));

	d->progressBar->setValue(0);

	d->checkBox_buildbifurcationmesh->setChecked(false);
	d->checkBox_loadlandmarks->setChecked(false);
	d->checkBox_loadcenterlines->setChecked(false);
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
		std::cout << "Will build bifurcation mesh" << std::endl;
		logic->WillBuildBifurcationMesh = true;
	}
	else
	{
		std::cout << "Will not build bifurcation mesh" << std::endl;
		logic->WillBuildBifurcationMesh = false;
	}

	logic->BuildCenterlinesMeshLogic();
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
			logic->DetectCenterlinesLogic(d->progressBar);
		}
		else
		{
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
			logic->BuildCenterlinesMeshLogic();
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