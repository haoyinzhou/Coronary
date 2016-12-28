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


class ORSliceStyle : public vtkInteractorStyleImage
{
public:
	static ORSliceStyle *New();
	vtkTypeMacro(ORSliceStyle, vtkInteractorStyleImage);

	ORSliceStyle()
	{
	//	this->widget = NULL;
		this->clModel = NULL;
		this->ObliqueReformat = NULL;
		this->obliqueImageSlicer = NULL;
		this->locator = vtkSmartPointer<vtkPointLocator>::New();

		this->pick = false;
		this->slide = false;
	//	clTube = NULL;
	}
	~ORSliceStyle()
	{
	}

	bool Pick(double picked[3])
	{
		int x = this->Interactor->GetEventPosition()[0];
		int y = this->Interactor->GetEventPosition()[1];

		this->FindPokedRenderer(x, y);
		if (this->CurrentRenderer == NULL) return false;

		vtkPropPicker *picker = vtkPropPicker::SafeDownCast(this->Interactor->GetPicker());
		if (picker == NULL) return false;

		// Pick at the mouse location provided by the interactor
		picker->Pick(x, y, 0.0, this->CurrentRenderer);

		// There could be other props assigned to this picker, so 
		// make sure we picked the image actor
		vtkAssemblyPath* path = picker->GetPath();
		bool validPick = false;

		if (path)
		{
			vtkCollectionSimpleIterator sit;
			path->InitTraversal(sit);
			vtkAssemblyNode *node;
			for (int i = 0; i < path->GetNumberOfItems() && !validPick; ++i)
			{
				node = path->GetNextNode(sit);
				if (obliqueImageSlicer == vtkImageSlice::SafeDownCast(node->GetViewProp()))
				{
					validPick = true;
				}
			}
		}

		if (!validPick)
			return false;

		// Get the world coordinates of the picked
		picker->GetPickPosition(picked);

		return true;
	}

	virtual void OnLeftButtonDown()
	{
	//	std::cout << "OnLeftButtonDown, cellid = " << cellid << std::endl;
	//	double pos[3] = {0.0, 0.0, 0.0};
	//	this->Pick(pos);

		if (!clModel || !ObliqueReformat || !obliqueImageSlicer)
			return;

		if (this->Pick(lastpickpos))
		{
			vtkPolyData *lumenPoly = vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(2));
			vtkPolyData *lumenCenter = vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(4));

			if (!lumenPoly || !lumenCenter) return;

			vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
			if (!clLumenRadius) return;

		//	std::cout << "lastpickpos = " << lastpickpos[0] << ", " << lastpickpos[1] << ", " << lastpickpos[2] << std::endl;

			locator = vtkSmartPointer<vtkPointLocator>::New();
			locator->SetDataSet(lumenPoly);
			locator->SetNumberOfPointsPerBucket(5);
			locator->AutomaticOn();
			locator->BuildLocator();
			double dist2;
			vtkIdType focalId = locator->FindClosestPointWithinRadius(4.0, lastpickpos, dist2);
		//	std::cout << "focalId = " << focalId << ", dist2 = " << dist2 << std::endl;

			if (focalId >= 0 && focalId < lumenPoly->GetNumberOfPoints())
			{
				vtkDoubleArray *paramArray = vtkDoubleArray::SafeDownCast(lumenPoly->GetPointData()->GetArray("Param"));
				if (!paramArray) return;

				double param[2];
				paramArray->GetTuple(focalId, param);

				focalParam[0] = vtkIdType((SmartCoronary::LongitudinalRefineSteps + 1) * param[0] + 0.5);
				focalParam[1] = vtkIdType(param[1] + 0.5) % clLumenRadius->GetNumberOfComponents();
				vtkIdType segmentId = this->ObliqueReformat->GetSegmentId();

				if (segmentId >= 0 && segmentId < clModel->GetNumberOfCells())
				{
					vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(segmentId, idlist);
					if (focalParam[0] >= 0 && focalParam[0] < idlist->GetNumberOfIds() &&
						focalParam[1] >= 0 && focalParam[1] < clLumenRadius->GetNumberOfComponents())
					{
						focalParam[0] = idlist->GetId(focalParam[0]);
						pick = true;
						if (!this->Interactor->GetControlKey())
							ObliqueReformat->UpdateImageOff();
					}
				}
			}
		}
	}

	virtual void OnLeftButtonUp()
	{
		if (pick)
		{
			//locator->Delete();
			pick = false;
			ObliqueReformat->UpdateImageOn();
		}
	}

	virtual void OnRightButtonDown()
	{
		if (this->ObliqueReformat)
		{
			vtkIdType segmentId = this->ObliqueReformat->GetSegmentId();
			vtkIdType pId = this->ObliqueReformat->GetPointId();
			vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
			clModel->GetCellPoints(segmentId, idlist);

			if (segmentId >= 0 && segmentId < clModel->GetNumberOfCells()
				&& pId >= 0 && pId < idlist->GetNumberOfIds())
			{
				slide = true;
				ObliqueReformat->UpdateImageOn();
			}
		}
		//this->Superclass::OnRightButtonDown();
	}

	virtual void OnRightButtonUp()
	{
		if (slide)
		{
			slide = false;
		}
		//this->Superclass::OnRightButtonUp();
	}
	

	virtual void OnMouseMove()
	{
		//std::cout << "slide = " << slide << ", pick = " << pick << std::endl;
		if (slide)
		{
			int eventpos[2], lasteventpos[2];
			this->Interactor->GetEventPosition(eventpos);
			this->Interactor->GetLastEventPosition(lasteventpos);
			int step = max(eventpos[0] - lasteventpos[0], eventpos[1] - lasteventpos[1]);

			if (step != 0)
			{
				ObliqueReformat->SetPointId(ObliqueReformat->GetPointId() + step);
				widget->GetRenderWindow()->Render();

		//		std::cout << "ObliqueReformat.pid = " << ObliqueReformat->GetPointId() << ", param = " << param[0] << ", " << param[1] << std::endl;
			}
		}

		if (pick)
		{
			if (this->Pick(pickpos))
			{
				double move[3];
				vtkMath::Subtract(pickpos, lastpickpos, move);
				
				int eventpos[2], lasteventpos[2];
				this->Interactor->GetEventPosition(eventpos);
				this->Interactor->GetLastEventPosition(lasteventpos);

				if (this->Interactor->GetShiftKey())
				{
					int step = eventpos[1] - lasteventpos[1];

					vtkDoubleArray *clRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Radius"));
					vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
					double newradius;
					newradius = clRadius->GetValue(focalParam[0]) * (1.0 + step * 0.02);
					if (newradius < 0.1) newradius = 0.1;
					else if (newradius > 10.0) newradius = 10.0;
					clRadius->SetValue(focalParam[0], newradius);
					for (int j = 0; j < clLumenRadius->GetNumberOfComponents(); j++)
					{
						newradius = clLumenRadius->GetComponent(focalParam[0], j) * (1.0 + step * 0.02);
						if (newradius < 0.1) newradius = 0.1;
						else if (newradius > 10.0) newradius = 10.0;
						clLumenRadius->SetComponent(focalParam[0], j, newradius);
					}
				}

				else
				{
					if (this->Interactor->GetControlKey())
					{
						double axis1[3], axis2[3];
						vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis1"));
						vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("Axis2"));
						clAxis1->GetTuple(focalParam[0], axis1);
						clAxis2->GetTuple(focalParam[0], axis2);
						double coord[3];
						clModel->GetPoint(focalParam[0], coord);
						for (int k = 0; k < 3; k++)
						{
							coord[k] -= (pickpos[0] - lastpickpos[0]) * axis1[k] + (pickpos[1] - lastpickpos[1]) * axis2[k];
						}
						clModel->GetPoints()->SetPoint(focalParam[0], coord);
					}
					else
					{						
						vtkDoubleArray *clArray;
						clArray = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
						vtkPolyData *lumenCenter = vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(4));
						double coord[3], dir[3];
						lumenCenter->GetPoint(0, coord);
						//std::cout << "lumenCenter coord = " << coord[0] << ", " << coord[1] << ", " << coord[2] << std::endl;
						vtkMath::Subtract(lastpickpos, coord, dir);
						vtkMath::Normalize(dir);
						double moveproj = vtkMath::Dot(move, dir);

					//	std::cout << "moving focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;

						double newradius = clArray->GetComponent(focalParam[0], focalParam[1]);
						newradius += moveproj;
						if (newradius < 0.1) newradius = 0.1;
						else if (newradius > 10.0) newradius = 10.0;
						clArray->SetComponent(focalParam[0], focalParam[1], newradius);

					//	std::cout << "moving2 focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;
					//	std::cout << clArray->GetComponent(focalParam[0], focalParam[1]) << std::endl;
					}
				}

				clModel->Modified();
				this->Interactor->Render();

				std::swap(lastpickpos, pickpos);
			}
		}

		this->Superclass::OnMouseMove();
	}

	virtual void OnKeyPress()
	{
		std::string key = this->Interactor->GetKeySym();

		if (key == "Up" || key == "Down")
		{
			int step = (key == "Up") ? 1 : -1;
			if (step != 0)
			{
				ObliqueReformat->SetPointId(ObliqueReformat->GetPointId() + step);
				widget->GetRenderWindow()->Render();
			}
			return;
		}

	//	this->Superclass::OnKeyPress();
		return;
	}

public:

	QVTKWidget* widget;

	vtkPolyData* clModel;
	ImageObliqueReformat* ObliqueReformat;
	vtkImageSlice *obliqueImageSlicer;	

	vtkSmartPointer<vtkPointLocator> locator;
	vtkIdType focalParam[2];
	double pickpos[3];
	double lastpickpos[3];

	bool pick;
	bool slide;
//	ExtendTubeFilter	* clTube;

};
vtkStandardNewMacro(ORSliceStyle);



QVesselEditingWidget::QVesselEditingWidget()
{
	widget1 = new QVTKWidget;
	widget2 = new QVTKWidget;

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addWidget(widget1, 1);
	layout->addWidget(widget2, 1);
	setLayout(layout);	

	ORSliceStyleCallback = vtkSmartPointer<ORSliceStyle>::New();
	ORSliceStyleCallback->widget = this->widget2;
	widget2->GetInteractor()->SetInteractorStyle(ORSliceStyleCallback);
}

QVesselEditingWidget::~QVesselEditingWidget()
{
	widget1->deleteLater();
	delete[] widget1;

	widget2->deleteLater();
	delete[] widget2;
}



void QVesselEditingWidget::setvisibleslot(bool f)
{
	std::cout << "set visible slot" << std::endl;
	this->setVisible(f);

	int parentHeight = this->height();

	widget1->setMinimumHeight(0.65 * parentHeight);
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
	std::cout << "SelectID = " << SelectID << std::endl;
	{
		vtkSmartPointer<vtkRendererCollection> rendercollection = widget1->GetRenderWindow()->GetRenderers();
		rendercollection->InitTraversal();
		for (vtkIdType i = 0; i < rendercollection->GetNumberOfItems(); i++)
		{
			vtkSmartPointer<vtkRenderer> thisrender = vtkRenderer::SafeDownCast(rendercollection->GetNextItem());
			widget1->GetRenderWindow()->RemoveRenderer(thisrender);
		}

		CurvedReformat = vtkSmartPointer<ImageCurvedReformat>::New();
		CurvedReformat->SetInputData(0, ImageData);
		CurvedReformat->SetInputData(1, clModel);
		CurvedReformat->SetSegmentId(SelectID);
		CurvedReformat->SetTwistIndex(0);
		CurvedReformat->Update();

		//vtkSmartPointer<vtkPolyData> poly;
		//SaveVTKImage(CurvedReformat->GetOutput(), "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output0.mha");
		//poly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(1));
		//SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output1.vtp");
		//poly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(2));
		//SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output2.vtp");
		//poly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(3));
		//SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output3.vtp");
		//poly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(4));
		//SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output4.vtp");
		//poly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(5));
		//SavePolyData(poly, "C:\\work\\Coronary_Slicer\\testdata\\CurvedReformat_output5.vtp");

		vtkSmartPointer<vtkRenderer> CurvedRenderer = vtkSmartPointer<vtkRenderer>::New();

		vtkSmartPointer<vtkImageSliceMapper> CurvedimageResliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
		CurvedimageResliceMapper->SetInputConnection(CurvedReformat->GetOutputPort(0));
		vtkSmartPointer<vtkImageSlice> CurvedimageSlice = vtkSmartPointer<vtkImageSlice>::New();
		CurvedimageSlice->SetMapper(CurvedimageResliceMapper);
		CurvedimageSlice->GetProperty()->SetColorWindow(1358);
		CurvedimageSlice->GetProperty()->SetColorLevel(-27);
		CurvedRenderer->AddViewProp(CurvedimageSlice);

		vtkSmartPointer<vtkPolyDataMapper> LumenContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		LumenContourMapper->SetInputConnection(CurvedReformat->GetOutputPort(2));
		vtkSmartPointer<vtkActor> LumenContourActor = vtkSmartPointer<vtkActor>::New();
		LumenContourActor->SetMapper(LumenContourMapper);
		LumenContourActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		LumenContourActor->GetProperty()->SetLineWidth(3.0f);
		LumenContourActor->GetProperty()->SetOpacity(0.6);
		LumenContourActor->PickableOff();
		CurvedRenderer->AddActor(LumenContourActor);

		vtkSmartPointer<vtkPolyDataMapper> clMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		clMapper->SetInputConnection(CurvedReformat->GetOutputPort(4));
		vtkSmartPointer<vtkActor> clActor = vtkSmartPointer<vtkActor>::New();
		clActor->SetMapper(clMapper);
		clActor->GetProperty()->SetColor(0.3, 0.4, 0.9);
		clActor->GetProperty()->SetLineWidth(3.0f);
		clActor->GetProperty()->SetOpacity(1.0);
		clActor->PickableOff();
		CurvedRenderer->AddActor(clActor);

		widget1->GetRenderWindow()->AddRenderer(CurvedRenderer);

		vtkCamera* camera = CurvedRenderer->GetActiveCamera();
		camera->ParallelProjectionOn();
		camera->SetPosition(0, 30, 1);
		camera->SetFocalPoint(0, 30, 0);
		camera->SetParallelScale(30);

		widget1->GetRenderWindow()->Render();
	}

	{
		vtkSmartPointer<vtkRendererCollection> rendercollection = widget2->GetRenderWindow()->GetRenderers();
		rendercollection->InitTraversal();
		for (vtkIdType i = 0; i < rendercollection->GetNumberOfItems(); i++)
		{
			vtkSmartPointer<vtkRenderer> thisrender = vtkRenderer::SafeDownCast(rendercollection->GetNextItem());
			widget2->GetRenderWindow()->RemoveRenderer(thisrender);
		}

		ObliqueReformat = vtkSmartPointer<ImageObliqueReformat>::New();
		ObliqueReformat->SetInputData(0, ImageData);
		ObliqueReformat->SetInputData(1, clModel);
		ObliqueReformat->SetSegmentId(SelectID);
		ObliqueReformat->SetPointId(0);
		ObliqueReformat->Update();

		vtkSmartPointer<vtkRenderer> ObliqueRenderer = vtkSmartPointer<vtkRenderer>::New();

		vtkSmartPointer<vtkImageSliceMapper> ObliqueimageResliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
		ObliqueimageResliceMapper->SetInputConnection(ObliqueReformat->GetOutputPort(0));
		vtkSmartPointer<vtkImageSlice> ObliqueimageSlice = vtkSmartPointer<vtkImageSlice>::New();
		ObliqueimageSlice->SetMapper(ObliqueimageResliceMapper);
		ObliqueimageSlice->GetProperty()->SetColorWindow(1358);
		ObliqueimageSlice->GetProperty()->SetColorLevel(-27);
		ObliqueRenderer->AddViewProp(ObliqueimageSlice);

		vtkSmartPointer<vtkPolyDataMapper> LumenContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		LumenContourMapper->SetInputConnection(ObliqueReformat->GetOutputPort(2));
		vtkSmartPointer<vtkActor> LumenContourActor = vtkSmartPointer<vtkActor>::New();
		LumenContourActor->SetMapper(LumenContourMapper);
		LumenContourActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		LumenContourActor->GetProperty()->SetLineWidth(3.0f);
		LumenContourActor->GetProperty()->SetOpacity(0.6);
		LumenContourActor->PickableOff();
		ObliqueRenderer->AddActor(LumenContourActor);

		vtkSmartPointer<vtkPolyDataMapper> clMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		clMapper->SetInputConnection(ObliqueReformat->GetOutputPort(4));
		vtkSmartPointer<vtkActor> clActor = vtkSmartPointer<vtkActor>::New();
		clActor->SetMapper(clMapper);
		clActor->GetProperty()->SetColor(0.3, 0.4, 0.9);
		clActor->GetProperty()->SetPointSize(6.0f);
		clActor->GetProperty()->SetOpacity(1.0);
		clActor->PickableOff();
		ObliqueRenderer->AddActor(clActor);

		widget2->GetRenderWindow()->AddRenderer(ObliqueRenderer);

		vtkCamera* camera = ObliqueRenderer->GetActiveCamera();
		camera->ParallelProjectionOn();
		camera->SetPosition(0, 30, 1);
		camera->SetFocalPoint(0, 30, 0);
		camera->SetParallelScale(30);

		widget2->GetRenderWindow()->Render();

		this->ORSliceStyleCallback->clModel = this->clModel;
		this->ORSliceStyleCallback->ObliqueReformat = this->ObliqueReformat;
		this->ORSliceStyleCallback->obliqueImageSlicer = ObliqueimageSlice;
	}



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

	std::cout << "force render slot done!" << std::endl;
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
	QDesktopWidget *desktop = QApplication::desktop();
	int screenWidth = desktop->width();
	int screenHeight = desktop->height();

	VesselEditingWidget->resize(screenWidth / 5, screenHeight / 1.2);
	VesselEditingWidget->move(0, 0);
	VesselEditingWidget->setWindowTitle("VesselEditingWidget");
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
		std::cout << "pickpixel = " << pickpixel[0] << ", " << pickpixel[1] << std::endl;

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

	if (logic->centerlineTube->GetOutput(0)->GetNumberOfCells() == 0)
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
	
	
	for (int i = 0; i < addedvesselpickobservertag.size(); i++)
		RenderWindowInteractorthreeD->RemoveObserver(addedvesselpickobservertag.at(i));
	addedvesselpickobservertag.clear();

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
	VesselPickCallBack->clmodel = logic->centerlineTube->GetOutput(0);
	VesselPickCallBack->lumenmodel = logic->centerlineTube->GetOutput(2);
	VesselPickCallBack->imagedata = logic->imageData;


	vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback> CtrlKeyPressedInteractionCallback = vtkSmartPointer<vtkCtrlKeyPressedInteractionCallback>::New();
	CtrlKeyPressedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyPressedInteractionCallback->VesselPicker = VesselPicker;
	CtrlKeyPressedInteractionCallback->VesselPickCallBack = VesselPickCallBack;
	CtrlKeyPressedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	CtrlKeyPressedInteractionCallback->clmodel = logic->centerlineTube->GetOutput(0);
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyPressEvent, CtrlKeyPressedInteractionCallback));

	vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback> CtrlKeyReleasedInteractionCallback = vtkSmartPointer<vtkCtrlKeyReleasedInteractionCallback>::New();
	CtrlKeyReleasedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	CtrlKeyReleasedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyReleaseEvent, CtrlKeyReleasedInteractionCallback));

	std::cout << "new ctrl observe has been set" << std::endl;

	return true;
}




bool qSlicerCoronaryMainModuleWidget::TestButtonFunc()
{
	std::cout << "TestButtonFunc begin" << std::endl;

	Q_D(qSlicerCoronaryMainModuleWidget);

	d->logic()->BuildCenterlinesMeshLogic();
	SetupKeyMouseObserver();

//	d->logic()->TestLogic();

	std::cout << "TestButtonFunc end" << std::endl;

	return true;
}
