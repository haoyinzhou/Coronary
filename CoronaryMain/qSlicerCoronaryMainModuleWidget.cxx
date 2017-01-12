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
		this->superwidget = NULL;
		this->widget = NULL;
		this->clModel = NULL;
		this->ObliqueReformat = NULL;
		this->obliqueImageSlicer = NULL;
		this->locator = vtkSmartPointer<vtkPointLocator>::New();

		this->pick = false;
		this->slide = false;
		this->pickdis = 0.0;

		this->pickedids[0] = 0;
		this->pickedids[1] = 0;

		for (int l = 0; l < 3; l++)
		{
			pickpos[l] = 0.0;
			lastpickpos[l] = 0.0;
		}
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

		vtkSmartPointer<vtkPropPicker> picker = vtkPropPicker::SafeDownCast(this->Interactor->GetPicker());
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

	bool GetNeighorClPoints(vector< vtkIdType >* ids_out, vector< double >* dis_out)
	{		
		ids_out->clear();
		dis_out->clear();

		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(pickedids[0], idlist);

		for (int i = (pickedids[1] - 1) >= 0 ? pickedids[1] - 1 : 0; i >= 0; i --)
		{
			if (idlist->GetId(pickedids[1]) == idlist->GetId(i))
				continue;

			double dis = abs(pickedids[1] - i);
			if (dis > superwidget->smoothclradius)
				break;

			ids_out->push_back(idlist->GetId(i));
			dis_out->push_back(dis);
		}	

		for (int i = (pickedids[1] + 1) < idlist->GetNumberOfIds() ? (pickedids[1] + 1) : idlist->GetNumberOfIds() - 1; i < idlist->GetNumberOfIds(); i ++)
		{
			if (idlist->GetId(pickedids[1]) == idlist->GetId(i))
				continue;

			double dis = abs(pickedids[1] - i);
			if (dis > superwidget->smoothclradius)
				break;

			ids_out->push_back(idlist->GetId(i));
			dis_out->push_back(dis);
		}

		return true;
	}

	virtual void OnLeftButtonDown()
	{
	//	std::cout << "OnLeftButtonDown, cellid = " << cellid << std::endl;
	//	double pos[3] = {0.0, 0.0, 0.0};
	//	this->Pick(pos);

		if (!clModel || !ObliqueReformat || !obliqueImageSlicer)
			return;

	//	clModel->BuildCells();
	//	vtkSmartPointer<vtkCell> thisline = clModel->GetCell(ObliqueReformat->GetSegmentId());
		//std::cout << "npts = " << thisline->GetNumberOfPoints() << std::endl;
	//	vtkSmartPointer<vtkIdList> idlistthisline = thisline->GetPointIds();
	//	std::cout << "npts = " << idlistthisline->GetNumberOfIds() << std::endl;
	//	std::cout << "pts[npts-1] = " << idlistthisline->GetId(idlistthisline->GetNumberOfIds() - 1) << std::endl;
	//	std::cout << "ObliqueReformatid = " << ObliqueReformat->GetPointId() << std::endl;

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

			vtkIdType focalId = locator->FindClosestPointWithinRadius(15.0, lastpickpos, this->pickdis);
		//	std::cout << "focalId = " << focalId << ", dist2 = " << dist2 << std::endl;

			if (focalId >= 0 && focalId < lumenPoly->GetNumberOfPoints())
			{
				vtkDoubleArray *paramArray = vtkDoubleArray::SafeDownCast(lumenPoly->GetPointData()->GetArray("Param"));
				if (!paramArray) return;
				double param[2];
				paramArray->GetTuple(focalId, param);

			//	std::cout << "param = " << param[0] << ", " << param[1] << std::endl;
				
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
						pickedids[0] = segmentId;
						pickedids[1] = focalParam[0];

					//	std::cout << "1 focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;

						// focalParam[0] = idlist->GetId(focalParam[0]);
						focalParam[0] = idlist->GetId(ObliqueReformat->GetPointId());

						pick = true;
						if (!this->Interactor->GetControlKey())
							ObliqueReformat->UpdateImageOff();

			//			std::cout << "2 focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;
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

				pickedids[0] = segmentId;
				pickedids[1] = pId;
			}
		}
		//this->Superclass::OnRightButtonDown();
	}

	virtual void OnRightButtonUp()
	{
		slide = false;
		//this->Superclass::OnRightButtonUp();
	}
	
	virtual void OnMouseMove()
	{
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
				{
					int	 imageDims[3];
					double imageOrigins[3];
					double imageSpacings[3];
					vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->CurvedReformat->GetOutput());
					reformatImage->GetDimensions(imageDims);
					reformatImage->GetOrigin(imageOrigins);
					reformatImage->GetSpacing(imageSpacings);
					double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					superwidget->CurvedReformatLine->SetPoint1(point1);
					superwidget->CurvedReformatLine->SetPoint2(point2);
					superwidget->widget1->GetRenderWindow()->Render();
				}
				{
					int	 imageDims[3];
					double imageOrigins[3];
					double imageSpacings[3];
					vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->stretchCurvedReformat->GetOutput());
					reformatImage->GetDimensions(imageDims);
					reformatImage->GetOrigin(imageOrigins);
					reformatImage->GetSpacing(imageSpacings);
					double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					superwidget->stretchCurvedReformatLine->SetPoint1(point1);
					superwidget->stretchCurvedReformatLine->SetPoint2(point2);
					superwidget->widget3->GetRenderWindow()->Render();
				}
				superwidget->send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(5)));
			}
		}

		else if (pick)
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

					vector<vtkIdType> neighorclid;
					vector<double> distance;
					GetNeighorClPoints(&neighorclid, &distance);

					for (int j = 0; j < clLumenRadius->GetNumberOfComponents(); j++)
					{
						newradius = clLumenRadius->GetComponent(focalParam[0], j);
						double moveproj = clLumenRadius->GetComponent(focalParam[0], j) * step * 0.02;
						newradius += moveproj;
						if (newradius < 0.1)
						{
							newradius = 0.1;
							moveproj = 0.0;
						}
						else if (newradius > 10.0)
						{
							newradius = 10.0;
							moveproj = 0.0;
						}
						clLumenRadius->SetComponent(focalParam[0], j, newradius);

						for (int i = 0; i < neighorclid.size(); i++)
						{
							double vd = clLumenRadius->GetComponent(neighorclid.at(i), j) - clLumenRadius->GetComponent(focalParam[0], j);
							if (vd * moveproj > 0.0) continue;

							double hd = distance.at(i) / superwidget->smoothclradius;
							double weight = exp(-hd * hd * 5);
							double weightedmove = weight * moveproj;
							double neighorradius = weightedmove + clLumenRadius->GetComponent(neighorclid.at(i), j);
							neighorradius = neighorradius < 0.1 ? 0.1 : neighorradius;
							neighorradius = neighorradius > 10.0 ? 10.0 : neighorradius;

							clLumenRadius->SetComponent(neighorclid.at(i), j, neighorradius);
						}

				//		surperwidget->send_lumenradiuschanged(focalParam[0], j, newradius);
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
						double coord[3], coordmove[3];
						clModel->GetPoint(focalParam[0], coord);
						for (int k = 0; k < 3; k++)
						{
							coordmove[k] = -(pickpos[0] - lastpickpos[0]) * axis1[k] - (pickpos[1] - lastpickpos[1]) * axis2[k];
							coord[k] += coordmove[k];
						}
						clModel->GetPoints()->SetPoint(focalParam[0], coord);

						superwidget->send_clcoordchanged(focalParam[0], coord[0], coord[1], coord[2]);
						{
							int	 imageDims[3];
							double imageOrigins[3];
							double imageSpacings[3];
							vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->CurvedReformat->GetOutput());
							reformatImage->GetDimensions(imageDims);
							reformatImage->GetOrigin(imageOrigins);
							reformatImage->GetSpacing(imageSpacings);
							double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
							double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
							superwidget->CurvedReformatLine->SetPoint1(point1);
							superwidget->CurvedReformatLine->SetPoint2(point2);
							superwidget->widget1->GetRenderWindow()->Render();
						}
						{
							int	 imageDims[3];
							double imageOrigins[3];
							double imageSpacings[3];
							vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->stretchCurvedReformat->GetOutput());
							reformatImage->GetDimensions(imageDims);
							reformatImage->GetOrigin(imageOrigins);
							reformatImage->GetSpacing(imageSpacings);
							double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
							double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
							superwidget->stretchCurvedReformatLine->SetPoint1(point1);
							superwidget->stretchCurvedReformatLine->SetPoint2(point2);
							superwidget->widget3->GetRenderWindow()->Render();
						}
					/*	vector<vtkIdType> neighorclid;
						vector<double> distance;
						GetNeighorClPoints(&neighorclid, &distance);
						for (int i = 0; i < neighorclid.size(); i ++)
						{
						//	double rd = distance.at(i) / superwidget->smoothclradius;
						//	double weight = 1.0 - rd * rd;
							double hd = distance.at(i) / superwidget->smoothclradius;
							double weight = exp(-hd * hd * 5);
							double wm[3];
							for (int k = 0; k < 3; k ++) wm[k] = weight * coordmove[k];
							clModel->GetPoint(neighorclid.at(i), coord);
							for (int k = 0; k < 3; k++) 	coord[k] += wm[k];
							clModel->GetPoints()->SetPoint(neighorclid.at(i), coord);
							superwidget->send_clcoordchanged(neighorclid.at(i), coord[0], coord[1], coord[2]);
						}
						*/
					}
					else
					{	
						if (this->pickdis < 4.0)
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
							if (newradius < 0.1)
							{
								newradius = 0.1;
								moveproj = 0.0;
							}
							else if (newradius > 10.0)
							{
								newradius = 10.0;
								moveproj = 0.0;
							}
							clArray->SetComponent(focalParam[0], focalParam[1], newradius);
							
							vector<vtkIdType> neighorclid;
							vector<double> distance;
							GetNeighorClPoints(&neighorclid, &distance);
	
							for (int i = 0; i < neighorclid.size(); i++)
							{
								double vd = clArray->GetComponent(neighorclid.at(i), focalParam[1]) - clArray->GetComponent(focalParam[0], focalParam[1]);
								if (vd * moveproj > 0.0) continue;

								double hd = distance.at(i) / superwidget->smoothclradius;
								double weight = exp(-hd * hd * 5);
								double weightedmove = weight * moveproj;
								double neighorradius = weightedmove + clArray->GetComponent(neighorclid.at(i), focalParam[1]);
								neighorradius = neighorradius < 0.1 ? 0.1 : neighorradius;
								neighorradius = neighorradius > 10.0 ? 10.0 : neighorradius;

								clArray->SetComponent(neighorclid.at(i), focalParam[1], neighorradius);
							}
							//	surperwidget->send_lumenradiuschanged(focalParam[0], focalParam[1], newradius);
						}
					}
				}

				clModel->Modified();
				this->Interactor->Render();
				superwidget->widget1->GetRenderWindow()->Render();
				superwidget->widget3->GetRenderWindow()->Render();

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
				this->Interactor->Render();
				{
					int	 imageDims[3];
					double imageOrigins[3];
					double imageSpacings[3];
					vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->CurvedReformat->GetOutput());
					reformatImage->GetDimensions(imageDims);
					reformatImage->GetOrigin(imageOrigins);
					reformatImage->GetSpacing(imageSpacings);
					double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					superwidget->CurvedReformatLine->SetPoint1(point1);
					superwidget->CurvedReformatLine->SetPoint2(point2);
					superwidget->widget1->GetRenderWindow()->Render();
				}
				{
					int	 imageDims[3];
					double imageOrigins[3];
					double imageSpacings[3];
					vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->stretchCurvedReformat->GetOutput());
					reformatImage->GetDimensions(imageDims);
					reformatImage->GetOrigin(imageOrigins);
					reformatImage->GetSpacing(imageSpacings);
					double point1[3] = { imageOrigins[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
					superwidget->stretchCurvedReformatLine->SetPoint1(point1);
					superwidget->stretchCurvedReformatLine->SetPoint2(point2);
					superwidget->widget3->GetRenderWindow()->Render();
				}
				superwidget->send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(5)));
			}
		}
		else if (key == "Left" || key == "Right")
		{
			int step = (key == "Right") ? 1 : -1;
			if (step != 0)
			{
				superwidget->CurvedReformat->SetTwistIndex(superwidget->CurvedReformat->GetTwistIndex() + step);
				superwidget->widget1->GetRenderWindow()->Render();
				superwidget->stretchCurvedReformat->SetTwistIndex(superwidget->CurvedReformat->GetTwistIndex());
				superwidget->widget3->GetRenderWindow()->Render();
			}
			superwidget->send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(superwidget->CurvedReformat->GetOutput(5)));
		}

		else if (key == "Control_L")
		{
			ObliqueReformat->UpdateImageOn();
		}
		else if (key == "g")
		{
			superwidget->send_detectlumen(ObliqueReformat->GetSegmentId());
		//	superwidget->send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(5)));
		//	superwidget->send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(superwidget->CurvedReformat->GetOutput(5)));
		}

	//	this->Superclass::OnKeyPress();
		return;
	}

	virtual void OnLeave()
	{
		pick = false;
		ObliqueReformat->UpdateImageOn();

		slide = false;
	}

public:
	QVesselEditingWidget* superwidget;
	QVTKWidget* widget;

	vtkPolyData* clModel;
	ImageObliqueReformat* ObliqueReformat;
	vtkImageSlice *obliqueImageSlicer;	

	vtkSmartPointer<vtkPointLocator> locator;
	vtkIdType focalParam[2];
	
	double pickpos[3];
	double lastpickpos[3];

	bool pick;
	double pickdis;
	bool slide;

	vtkIdType pickedids[2]; // pickedids[0] segment id; pickedids[1] point id in this segment (not the real pid)
};
vtkStandardNewMacro(ORSliceStyle);


class CRRotateStyle : public vtkInteractorStyleImage
{
public:
	static CRRotateStyle *New();
	vtkTypeMacro(CRRotateStyle, vtkInteractorStyleImage);
	
	CRRotateStyle()
	{
		this->superwidget = NULL;
		this->widget = NULL;
		this->clModel = NULL;
		this->CurvedReformat = NULL;
		this->curvedImageSlicer = NULL;
		this->locator = vtkSmartPointer<vtkPointLocator>::New();
		this->focalId = -1;

		this->pick = false;
		this->rotate = false;
		this->pickdis = 0;

		this->pickedids[0] = 0;
		this->pickedids[1] = 0;

		for (int l = 0; l < 3; l++)
		{
			pickpos[l] = 0.0;
			lastpickpos[l] = 0.0;
		}
	}

	~CRRotateStyle()
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
				if (curvedImageSlicer == vtkImageSlice::SafeDownCast(node->GetViewProp()))
					validPick = true;
			}
		}
		if (!validPick)
			return false;
		// Get the world coordinates of the picked
		picker->GetPickPosition(picked);
		
		return true;
	}


	bool GetNeighorClPoints(vector< vtkIdType >* ids_out, vector< double >* dis_out)
	{
		ids_out->clear();
		dis_out->clear();

		vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
		clModel->GetCellPoints(pickedids[0], idlist);

		for (int i = (pickedids[1] - 1) >= 0 ? pickedids[1] - 1 : 0; i >= 0; i--)
		{
			if (idlist->GetId(pickedids[1]) == idlist->GetId(i))
				continue;

			double dis = abs(pickedids[1] - i);
			if (dis > superwidget->smoothclradius)
				break;

			ids_out->push_back(idlist->GetId(i));
			dis_out->push_back(dis);
		}

		for (int i = (pickedids[1] + 1) < idlist->GetNumberOfIds() ? (pickedids[1] + 1) : idlist->GetNumberOfIds() - 1; i < idlist->GetNumberOfIds(); i++)
		{
			if (idlist->GetId(pickedids[1]) == idlist->GetId(i))
				continue;

			double dis = abs(pickedids[1] - i);
			if (dis > superwidget->smoothclradius)
				break;

			ids_out->push_back(idlist->GetId(i));
			dis_out->push_back(dis);
		}

		return true;
	}
	
	virtual void OnLeftButtonDown()
	{
		if (!clModel || !CurvedReformat || !curvedImageSlicer)
			return;

		if (this->Pick(lastpickpos))
		{
			vtkPolyData *lumenPoly = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(2));
			vtkPolyData *lumenCenter = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(4));

			if (!lumenPoly || !lumenCenter) return;

			vtkDoubleArray *clLumenRadius = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
			if (!clLumenRadius) return;

			locator = vtkSmartPointer<vtkPointLocator>::New();
			locator->SetDataSet(lumenPoly);
			locator->SetNumberOfPointsPerBucket(5);
			locator->AutomaticOn();
			locator->BuildLocator();

			this->focalId = locator->FindClosestPointWithinRadius(4.0, lastpickpos, this->pickdis);
			//	std::cout << "focalId = " << focalId << ", dist2 = " << dist2 << std::endl;

			if (this->focalId >= 0 && this->focalId < lumenPoly->GetNumberOfPoints())
			{
				vtkDoubleArray *paramArray = vtkDoubleArray::SafeDownCast(lumenPoly->GetPointData()->GetArray("Param"));
				if (!paramArray) return;

				double param[2];
				paramArray->GetTuple(focalId, param);

				focalParam[0] = vtkIdType((SmartCoronary::LongitudinalRefineSteps + 1) * param[0] + 0.5);
				focalParam[1] = vtkIdType(param[1] + 0.5) % clLumenRadius->GetNumberOfComponents();
				vtkIdType segmentId = this->CurvedReformat->GetSegmentId();

				if (segmentId >= 0 && segmentId < clModel->GetNumberOfCells())
				{
					vtkSmartPointer<vtkIdList> idlist = vtkSmartPointer<vtkIdList>::New();
					clModel->GetCellPoints(segmentId, idlist);
					if (focalParam[0] >= 0 && focalParam[0] < idlist->GetNumberOfIds() &&
						focalParam[1] >= 0 && focalParam[1] < clLumenRadius->GetNumberOfComponents())
					{
						pickedids[0] = segmentId;
						pickedids[1] = focalParam[0];

						focalParam[0] = idlist->GetId(focalParam[0]);
						pick = true;
					//	if (!this->Interactor->GetControlKey())
					//		CurvedReformat->UpdateImageOff();

				//		std::cout << "this->focalId = " << this->focalId << "focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;
					}
				}
			}
		}
	}
	
	virtual void OnLeftButtonUp()
	{
		if (pick)
		{
			pick = false;
			this->focalId = -1;
		//	CurvedReformat->UpdateImageOn();
		}
	}
	
	virtual void OnRightButtonDown()
	{
		if (CurvedReformat && CurvedReformat->GetSegmentId() >= 0)
		{
			rotate = true;
			CurvedReformat->UpdateImageOn();
		}
		return;
		//this->Superclass::OnRightButtonDown();
	}
	virtual void OnRightButtonUp()
	{
		rotate = false;
		//this->Superclass::OnRightButtonUp();
	}

	virtual void OnMouseMove()
	{
		if (rotate)
		{
			int eventpos[2], lasteventpos[2];
			this->Interactor->GetEventPosition(eventpos);
			this->Interactor->GetLastEventPosition(lasteventpos);
			int step = max(eventpos[0] - lasteventpos[0], eventpos[1] - lasteventpos[1]);

			if (step != 0)
			{
				CurvedReformat->SetTwistIndex(CurvedReformat->GetTwistIndex() + step);
				this->Interactor->Render();

				superwidget->stretchCurvedReformat->SetTwistIndex(CurvedReformat->GetTwistIndex());
				superwidget->widget3->GetRenderWindow()->Render();

				superwidget->send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(5)));
			}
		}
		else if (pick)
		{
			if (Pick(pickpos))
			{
			//	this->clModel->Modified();

				double move[3];
				vtkMath::Subtract(pickpos, lastpickpos, move);

				int eventpos[2], lasteventpos[2];
				this->Interactor->GetEventPosition(eventpos);
				this->Interactor->GetLastEventPosition(lasteventpos);
		//		std::cout << "eventpos = " << eventpos[0] << ", " << eventpos[1] << std::endl;
		//		std::cout << "lasteventpos = " << lasteventpos[0] << ", " << lasteventpos[1] << std::endl;
				
				if (this->Interactor->GetShiftKey())
				{
				}
				else
				{
					double coord[3], dir[3];
					vtkPolyData *lumenCenter = vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(4));
					lumenCenter->GetPoint(focalId / 2, coord);
					vtkMath::Subtract(lastpickpos, coord, dir);
					vtkMath::Normalize(dir);
					double moveproj = vtkMath::Dot(move, dir);

					if (this->pickdis < 4.0)
					{
						vtkDoubleArray *clArray = vtkDoubleArray::SafeDownCast(clModel->GetPointData()->GetArray("LumenRadius"));
					//	std::cout << "moving focalParam = " << focalParam[0] << ", " << focalParam[1] << std::endl;

						double newradius = clArray->GetComponent(focalParam[0], focalParam[1]);
						newradius += moveproj;
						if (newradius < 0.1)
						{
							newradius = 0.1;
							moveproj = 0.0;
						}
						else if (newradius > 10.0)
						{
							newradius = 10.0;
							moveproj = 0.0;
						}
						clArray->SetComponent(focalParam[0], focalParam[1], newradius);

						vector<vtkIdType> neighorclid;
						vector<double> distance;
						GetNeighorClPoints(&neighorclid, &distance);

						for (int i = 0; i < neighorclid.size(); i++)
						{
							double vd = clArray->GetComponent(neighorclid.at(i), focalParam[1]) - clArray->GetComponent(focalParam[0], focalParam[1]);
							if (vd * moveproj > 0.0) continue;

							double hd = distance.at(i) / superwidget->smoothclradius;
							double weight = exp(-hd * hd * 5);
							double weightedmove = weight * moveproj;
							double neighorradius = weightedmove + clArray->GetComponent(neighorclid.at(i), focalParam[1]);
							neighorradius = neighorradius < 0.1 ? 0.1 : neighorradius;
							neighorradius = neighorradius > 10.0 ? 10.0 : neighorradius;

							clArray->SetComponent(neighorclid.at(i), focalParam[1], neighorradius);
						}
						
							//	surperwidget->send_lumenradiuschanged(focalParam[0], focalParam[1], newradius);
				//		SmoothLumenRadius(clModel, this->CurvedReformat->GetSegmentId(), 2);
					}				
				}			
				this->clModel->Modified();
				this->Interactor->Render();
				superwidget->widget2->GetRenderWindow()->Render();
				superwidget->widget3->GetRenderWindow()->Render();

				std::swap(lastpickpos, pickpos);	
			}
		}		
		
		this->Superclass::OnMouseMove();
	}
	
	virtual void OnKeyPress()
	{
		std::string key = this->Interactor->GetKeySym();

		if (key == "Left" || key == "Right")
		{
			int step = (key == "Right") ? 1 : -1;
			if (step != 0)
			{
				CurvedReformat->SetTwistIndex(CurvedReformat->GetTwistIndex() + step);
				this->Interactor->Render();
				superwidget->stretchCurvedReformat->SetTwistIndex(CurvedReformat->GetTwistIndex());
				superwidget->widget3->GetRenderWindow()->Render();

				superwidget->send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(5)));

			}
		}
		else if (key == "Up" || key == "Down")
		{
			int step = (key == "Up") ? 1 : -1;
			if (step != 0)
			{
				superwidget->ObliqueReformat->SetPointId(superwidget->ObliqueReformat->GetPointId() + step);
				this->Interactor->Render();

				int	 imageDims[3];
				double imageOrigins[3];
				double imageSpacings[3];
				vtkImageData* reformatImage = vtkImageData::SafeDownCast(superwidget->CurvedReformat->GetOutput());
				reformatImage->GetDimensions(imageDims);
				reformatImage->GetOrigin(imageOrigins);
				reformatImage->GetSpacing(imageSpacings);
				double point1[3] = { imageOrigins[0], imageOrigins[1] + superwidget->ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
				double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + superwidget->ObliqueReformat->GetPointId() * imageSpacings[1], imageOrigins[2] };
				superwidget->CurvedReformatLine->SetPoint1(point1);
				superwidget->CurvedReformatLine->SetPoint2(point2);

				this->Interactor->Render();
				superwidget->widget2->GetRenderWindow()->Render();
				superwidget->widget3->GetRenderWindow()->Render();

				superwidget->send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(superwidget->ObliqueReformat->GetOutput(5)));
			}
		}
		else if (key == "g")
		{
			superwidget->send_detectlumen(CurvedReformat->GetSegmentId());
		//	superwidget->send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(5)));
		//	superwidget->send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(superwidget->ObliqueReformat->GetOutput(5)));
		}


		//	this->Superclass::OnKeyPress();
		return;
	}

	virtual void OnLeave()
	{
		pick = false;
		rotate = false;
	}

public:

	QVesselEditingWidget* superwidget;
	QVTKWidget* widget;

	vtkPolyData* clModel;
	ImageCurvedReformat* CurvedReformat;
	vtkImageSlice *curvedImageSlicer;

	vtkSmartPointer<vtkPointLocator> locator;
	vtkIdType focalParam[2];
	vtkIdType focalId;
	
	double pickpos[3];
	double lastpickpos[3];
	double pickdis;

	bool pick;
	bool rotate;

	vtkIdType pickedids[2]; // pickedids[0] segment id; pickedids[1] point id in this segment (not the real pid)

};
vtkStandardNewMacro(CRRotateStyle);


class SCRRotateStyle : public vtkInteractorStyleImage
{
public:
	static SCRRotateStyle *New();
	vtkTypeMacro(SCRRotateStyle, vtkInteractorStyleImage);

	SCRRotateStyle()
	{
		this->superwidget = NULL;
		this->widget = NULL;
	}
	~SCRRotateStyle(){}


public:
	QVesselEditingWidget* superwidget;
	QVTKWidget* widget;

};
vtkStandardNewMacro(SCRRotateStyle);


QVesselEditingWidget::QVesselEditingWidget()
{
	this->widget1 = new QVTKWidget;
	this->widget2 = new QVTKWidget;
	this->widget3 = new QVTKWidget;

	QHBoxLayout* layout1 = new QHBoxLayout;
	layout1->addWidget(widget1);
	layout1->addWidget(widget3);

	QVBoxLayout *layout = new QVBoxLayout;
	layout->addLayout(layout1);
	layout->addWidget(widget2, 1);
	this->setLayout(layout);
		
	this->ORSliceStyleCallback = vtkSmartPointer<ORSliceStyle>::New();
	this->ORSliceStyleCallback->superwidget = this;
	this->ORSliceStyleCallback->widget = this->widget2;
	this->widget2->GetInteractor()->SetInteractorStyle(ORSliceStyleCallback);

	this->CRRotateStyleCallback = vtkSmartPointer<CRRotateStyle>::New();
	this->CRRotateStyleCallback->superwidget = this;
	this->CRRotateStyleCallback->widget = this->widget1;
	this->widget1->GetInteractor()->SetInteractorStyle(CRRotateStyleCallback);

	this->Visiblity = false;
	this->smoothclradius = 5.5;

	this->CurvedReformatLine = vtkSmartPointer<vtkLineSource>::New();
	this->stretchCurvedReformatLine = vtkSmartPointer<vtkLineSource>::New();
}

QVesselEditingWidget::~QVesselEditingWidget()
{
	widget1->deleteLater();
	delete[] widget1;

	widget2->deleteLater();
	delete[] widget2;

	widget3->deleteLater();
	delete[] widget3;

}

void QVesselEditingWidget::closeEvent(QCloseEvent *event)
{
	this->Visiblity = false;
	emit widgetclosedsignal();
	QVTKWidget::closeEvent(event);
}

void QVesselEditingWidget::setvisibleslot(bool f)
{
	this->setVisible(f);

	int parentHeight = this->height();
	widget1->setMinimumHeight(0.5 * parentHeight);
	widget3->setMinimumHeight(0.5 * parentHeight);
	widget2->setMinimumHeight(0.35 * parentHeight);

	if (Visiblity == false && f == true)
	{
		emit removemouseobserveratmainwidget();
	}

	Visiblity = f;
}


void QVesselEditingWidget::setselectidslot(vtkIdType id)
{
	this->SelectID = id;
}

void QVesselEditingWidget::setclmodelslot(vtkPolyData* cl)
{
	this->clModel = cl;
}

void QVesselEditingWidget::setimagedataslot(vtkImageData* p)
{
	this->ImageData = p;
}

void QVesselEditingWidget::resetslot()
{
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
		CurvedReformat->UpdateImageOn();
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
		
		int	 imageDims[3];
		double imageOrigins[3];
		double imageSpacings[3];
		vtkImageData* reformatImage = vtkImageData::SafeDownCast(CurvedReformat->GetOutput());
		reformatImage->GetDimensions(imageDims);
		reformatImage->GetOrigin(imageOrigins);
		reformatImage->GetSpacing(imageSpacings);
		double point1[3] = { imageOrigins[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		CurvedReformatLine->SetPoint1(point1);
		CurvedReformatLine->SetPoint2(point2);

		vtkSmartPointer<vtkPolyDataMapper> curveReformatLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		curveReformatLineMapper->SetInputConnection(this->CurvedReformatLine->GetOutputPort());
		vtkSmartPointer<vtkActor> curveReformatLineActor = vtkSmartPointer<vtkActor>::New();
		curveReformatLineActor->SetMapper(curveReformatLineMapper);
		curveReformatLineActor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		curveReformatLineActor->GetProperty()->SetLineWidth(3.0f);
		curveReformatLineActor->GetProperty()->SetOpacity(0.6);
		curveReformatLineActor->PickableOff();
		CurvedRenderer->AddActor(curveReformatLineActor);
		
		widget1->GetRenderWindow()->AddRenderer(CurvedRenderer);

		vtkCamera* camera = CurvedRenderer->GetActiveCamera();
		camera->ParallelProjectionOn();
		camera->SetPosition(0, 30, 1);
		camera->SetFocalPoint(0, 30, 0);
		camera->SetParallelScale(30);
	
		widget1->GetRenderWindow()->Render();

		this->CRRotateStyleCallback->clModel = this->clModel;
		this->CRRotateStyleCallback->CurvedReformat = this->CurvedReformat;
		this->CRRotateStyleCallback->curvedImageSlicer = CurvedimageSlice;

		send_updatecurvedslicesignal(vtkPolyData::SafeDownCast(CurvedReformat->GetOutput(5)));

//		CurvedReformat->GetOutputPort(5);
	}

	{
		vtkSmartPointer<vtkRendererCollection> rendercollection = widget3->GetRenderWindow()->GetRenderers();
		rendercollection->InitTraversal();
		for (vtkIdType i = 0; i < rendercollection->GetNumberOfItems(); i++)
		{
			vtkSmartPointer<vtkRenderer> thisrender = vtkRenderer::SafeDownCast(rendercollection->GetNextItem());
			widget3->GetRenderWindow()->RemoveRenderer(thisrender);
		}

		stretchCurvedReformat = vtkSmartPointer<ImageStretchCurvedReformat>::New();
		stretchCurvedReformat->SetInputData(0, ImageData);
		stretchCurvedReformat->SetInputData(1, clModel);
		stretchCurvedReformat->SetSegmentId(SelectID);
		stretchCurvedReformat->SetTwistIndex(0);
		stretchCurvedReformat->UpdateImageOn();
		stretchCurvedReformat->Update();

		vtkSmartPointer<vtkRenderer> stretchCurvedRenderer = vtkSmartPointer<vtkRenderer>::New();

		vtkSmartPointer<vtkImageSliceMapper> stretchCurvedimageResliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
		stretchCurvedimageResliceMapper->SetInputConnection(stretchCurvedReformat->GetOutputPort(0));
		vtkSmartPointer<vtkImageSlice> stretchCurvedimageSlice = vtkSmartPointer<vtkImageSlice>::New();
		stretchCurvedimageSlice->SetMapper(stretchCurvedimageResliceMapper);
		stretchCurvedimageSlice->GetProperty()->SetColorWindow(1358);
		stretchCurvedimageSlice->GetProperty()->SetColorLevel(-27);
		stretchCurvedRenderer->AddViewProp(stretchCurvedimageSlice);

		vtkSmartPointer<vtkPolyDataMapper> LumenContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		LumenContourMapper->SetInputConnection(stretchCurvedReformat->GetOutputPort(2));
		vtkSmartPointer<vtkActor> LumenContourActor = vtkSmartPointer<vtkActor>::New();
		LumenContourActor->SetMapper(LumenContourMapper);
		LumenContourActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		LumenContourActor->GetProperty()->SetLineWidth(3.0f);
		LumenContourActor->GetProperty()->SetOpacity(0.6);
		LumenContourActor->PickableOff();
		stretchCurvedRenderer->AddActor(LumenContourActor);

		vtkSmartPointer<vtkPolyDataMapper> clMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		clMapper->SetInputConnection(stretchCurvedReformat->GetOutputPort(4));
		vtkSmartPointer<vtkActor> clActor = vtkSmartPointer<vtkActor>::New();
		clActor->SetMapper(clMapper);
		clActor->GetProperty()->SetColor(0.3, 0.4, 0.9);
		clActor->GetProperty()->SetLineWidth(3.0f);
		clActor->GetProperty()->SetOpacity(1.0);
		clActor->PickableOff();
		stretchCurvedRenderer->AddActor(clActor);

/*		int	 imageDims[3];
		double imageOrigins[3];
		double imageSpacings[3];
		vtkImageData* reformatImage = vtkImageData::SafeDownCast(CurvedReformat->GetOutput());
		reformatImage->GetDimensions(imageDims);
		reformatImage->GetOrigin(imageOrigins);
		reformatImage->GetSpacing(imageSpacings);
		double point1[3] = { imageOrigins[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		CurvedReformatLine->SetPoint1(point1);
		CurvedReformatLine->SetPoint2(point2);

		vtkSmartPointer<vtkPolyDataMapper> curveReformatLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		curveReformatLineMapper->SetInputConnection(this->CurvedReformatLine->GetOutputPort());
		vtkSmartPointer<vtkActor> curveReformatLineActor = vtkSmartPointer<vtkActor>::New();
		curveReformatLineActor->SetMapper(curveReformatLineMapper);
		curveReformatLineActor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		curveReformatLineActor->GetProperty()->SetLineWidth(3.0f);
		curveReformatLineActor->GetProperty()->SetOpacity(0.6);
		curveReformatLineActor->PickableOff();
		CurvedRenderer->AddActor(curveReformatLineActor);
*/
		
		int	 imageDims[3];
		double imageOrigins[3];
		double imageSpacings[3];
		vtkImageData* reformatImage = vtkImageData::SafeDownCast(stretchCurvedReformat->GetOutput());
		reformatImage->GetDimensions(imageDims);
		reformatImage->GetOrigin(imageOrigins);
		reformatImage->GetSpacing(imageSpacings);
		double point1[3] = { imageOrigins[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		double point2[3] = { imageOrigins[0] + (imageDims[0] - 1)*imageSpacings[0], imageOrigins[1] + 0 * imageSpacings[1], imageOrigins[2] };
		stretchCurvedReformatLine->SetPoint1(point1);
		stretchCurvedReformatLine->SetPoint2(point2);

		vtkSmartPointer<vtkPolyDataMapper> stretchcurveReformatLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		stretchcurveReformatLineMapper->SetInputConnection(this->stretchCurvedReformatLine->GetOutputPort());
		vtkSmartPointer<vtkActor> stretchcurveReformatLineActor = vtkSmartPointer<vtkActor>::New();
		stretchcurveReformatLineActor->SetMapper(stretchcurveReformatLineMapper);
		stretchcurveReformatLineActor->GetProperty()->SetColor(0.0, 0.0, 1.0);
		stretchcurveReformatLineActor->GetProperty()->SetLineWidth(3.0f);
		stretchcurveReformatLineActor->GetProperty()->SetOpacity(0.6);
		stretchcurveReformatLineActor->PickableOff();
		stretchCurvedRenderer->AddActor(stretchcurveReformatLineActor);
		
		widget3->GetRenderWindow()->AddRenderer(stretchCurvedRenderer);

		vtkCamera* camera = stretchCurvedRenderer->GetActiveCamera();
		camera->ParallelProjectionOn();
		camera->SetPosition(0, 30, 1);
		camera->SetFocalPoint(0, 30, 0);
		camera->SetParallelScale(30);

		widget3->GetRenderWindow()->Render();
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

		send_updateobliqueslicesignal(vtkPolyData::SafeDownCast(ObliqueReformat->GetOutput(5)));
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

}

void QVesselEditingWidget::simplerenderslot()
{
	this->widget1->GetRenderWindow()->Render();
	this->widget2->GetRenderWindow()->Render();
	this->widget3->GetRenderWindow()->Render();
}


void QVesselEditingWidget::send_clcoordchanged(vtkIdType pointid, double x, double y, double z)
{
	emit clcoordchanged(pointid, x, y, z);
}
void QVesselEditingWidget::send_lumenradiuschanged(vtkIdType pointid, vtkIdType lumenpointid, double radius)
{
	emit lumenradiuschanged(pointid, lumenpointid, radius);
}

void QVesselEditingWidget::send_detectlumen(vtkIdType sid)
{
	emit detectlumensignal(sid);
}

void QVesselEditingWidget::send_updateobliqueslicesignal(vtkPolyData* p)
{
	emit updateobliqueslicesignal(p);
}

void QVesselEditingWidget::send_updatecurvedslicesignal(vtkPolyData* p)
{
	emit updatecurvedslicesignal(p);
}



void qSlicerCoronaryMainModuleWidget::send_visibilitychanged(bool f)
{
	emit visibilitychanged(f);
}
void qSlicerCoronaryMainModuleWidget::send_selectidchanged(vtkIdType id)
{
	emit selectidchanged(id);
}
void qSlicerCoronaryMainModuleWidget::send_clmodelchanged(vtkPolyData* cl)
{
	emit clmodelchanged(cl);
}
void qSlicerCoronaryMainModuleWidget::send_imagedatachanged(vtkImageData* p)
{
	emit imagedatachanged(p);
}
void qSlicerCoronaryMainModuleWidget::send_resetsignal()
{
	emit resetsignal();
}
void qSlicerCoronaryMainModuleWidget::send_forcerendersignal()
{
	emit forcerendersignal();
}
void qSlicerCoronaryMainModuleWidget::send_simplerendersignal()
{
	emit simplerendersignal();
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
	this->SelectedVesselID = -1;
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

	ObliqueSlicePolydata = vtkSmartPointer<vtkPolyData>::New();
	CurvedSlicePolydata = vtkSmartPointer<vtkPolyData>::New();

	VesselEditingWidget = new QVesselEditingWidget;
	QDesktopWidget *desktop = QApplication::desktop();
	int screenWidth = desktop->width();
	int screenHeight = desktop->height();

	VesselEditingWidget->resize(screenWidth / 5, screenHeight / 1.2);
	VesselEditingWidget->move(0, 0);
	VesselEditingWidget->setWindowTitle("Vessel Editing Widget");
	VesselEditingWidget->setVisible(false);

	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
	this->threeDView = layoutManager->threeDWidget(0)->threeDView()->VTKWidget();
	this->RenderWindowInteractorthreeD = threeDView->GetInteractor();
	
	
	connect(d->DetectLandmarks, SIGNAL(clicked()), this, SLOT(DetectLandmarksButtonFunc()));
	connect(d->DetectCenterlines, SIGNAL(clicked()), this, SLOT(DetectCenterlinesButtonFunc()));
	connect(d->DetectLumen, SIGNAL(clicked()), this, SLOT(DetectLumenButtonFunc()));
	connect(d->SaveLandmarks, SIGNAL(clicked()), this, SLOT(SaveLandmarksButtonFunc()));
	connect(d->SaveCenterlines, SIGNAL(clicked()), this, SLOT(SaveCenterlinesButtonFunc()));
	connect(d->MRMLNodeReadVolumn, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(SetVolumn(vtkMRMLNode*)));
	connect(d->checkBox_buildbifurcationmesh, SIGNAL(stateChanged(int)), this, SLOT(SetCheckBoxBuildBifurcationMesh(int)));

	connect(this, SIGNAL(visibilitychanged(bool)), VesselEditingWidget, SLOT(setvisibleslot(bool)));
	connect(this, SIGNAL(selectidchanged(vtkIdType)), VesselEditingWidget, SLOT(setselectidslot(vtkIdType)));
	connect(this, SIGNAL(clmodelchanged(vtkPolyData*)), VesselEditingWidget, SLOT(setclmodelslot(vtkPolyData*)));
	connect(this, SIGNAL(imagedatachanged(vtkImageData*)), VesselEditingWidget, SLOT(setimagedataslot(vtkImageData*)));
	connect(this, SIGNAL(resetsignal(void)), VesselEditingWidget, SLOT(resetslot()));
	connect(this, SIGNAL(forcerendersignal(void)), VesselEditingWidget, SLOT(forcerenderslot()));
	connect(this, SIGNAL(simplerendersignal()), VesselEditingWidget, SLOT(simplerenderslot()));
	connect(VesselEditingWidget, SIGNAL(clcoordchanged(vtkIdType, double, double, double)), this, SLOT(setclcoordslot(vtkIdType, double, double, double)));
	connect(VesselEditingWidget, SIGNAL(lumenradiuschanged(vtkIdType, vtkIdType, double)), this, SLOT(setlumenradiusslot(vtkIdType, vtkIdType, double)));
	connect(VesselEditingWidget, SIGNAL(removemouseobserveratmainwidget()), this, SLOT(removemouseobserverslot()));
	connect(VesselEditingWidget, SIGNAL(detectlumensignal(vtkIdType)), this, SLOT(detectlumenslot(vtkIdType)));
	connect(VesselEditingWidget, SIGNAL(widgetclosedsignal()), this, SLOT(vesseleditingwidgetclosedslot()));
	connect(VesselEditingWidget, SIGNAL(updateobliqueslicesignal(vtkPolyData*)), this, SLOT(updateobliquesliceslot(vtkPolyData*)));
	connect(VesselEditingWidget, SIGNAL(updatecurvedslicesignal(vtkPolyData*)), this, SLOT(updatecurvedsliceslot(vtkPolyData*)));

	connect(d->pushButtonTest, SIGNAL(clicked()), this, SLOT(TestButtonFunc()));

	d->progressBar->setValue(0);
	d->checkBox_buildbifurcationmesh->setChecked(false);
	d->checkBox_loadlandmarks->setChecked(false);
	d->checkBox_loadcenterlines->setChecked(false);

	addedselectedclnode.clear();
	addedobliqueandcurvednode.clear();
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
			for (int i = SmartCoronary::LEFT_CORONARY_OSTIUM; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
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
			for (int i = SmartCoronary::LEFT_CORONARY_OSTIUM; i < SmartCoronary::NUMBER_OF_LVCOR_LANDMARKS; i++)
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
			this->send_visibilitychanged(false);

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
			logic->centerlineModel_display = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel = vtkSmartPointer<vtkPolyData>::New();
			logic->LumenModel_display = vtkSmartPointer<vtkPolyData>::New();

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

		this->send_visibilitychanged(false);

		d->progressBar->setValue(100);
	}
	return true;
}

bool qSlicerCoronaryMainModuleWidget::DetectLumenButtonFunc()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	this->send_visibilitychanged(false);

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

	this->SelectedVesselID = -1;

	if (addedselectedclnode.size() != 0)
	{
		for (int i = 0; i < addedselectedclnode.size(); i++)
			logic->GetMRMLScene()->RemoveNode(addedselectedclnode.at(i));
		addedselectedclnode.clear();
	}
	return true;
}

bool qSlicerCoronaryMainModuleWidget
::RemoveAllObliqueandCurvedSlice()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (addedobliqueandcurvednode.size() != 0)
	{
		for (int i = 0; i < addedobliqueandcurvednode.size(); i++)
			logic->GetMRMLScene()->RemoveNode(addedobliqueandcurvednode.at(i));
		addedobliqueandcurvednode.clear();
	}
	return true;
}



bool qSlicerCoronaryMainModuleWidget
::ShowSelectedVesselThreeD(vtkIdType cellid)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();
	
	vtkSmartPointer<vtkIdTypeArray> centerlineSelectId = vtkSmartPointer<vtkIdTypeArray>::New();
	centerlineSelectId->SetName("SegmentId");
	centerlineSelectId->SetNumberOfValues(1);
	centerlineSelectId->SetValue(0, cellid);

	if (cellid >= logic->centerlineModel->GetNumberOfCells())
		cellid = 0;

	RemoveAllSelectedVesselThreeD();
	this->SelectedVesselID = cellid;

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
	{
		vtkMRMLNode* thisaddednode;

		SelectedClNode = vtkSmartPointer< vtkMRMLModelNode >::New();
		SelectedClDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();

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
	}
	
	return true;
}

bool qSlicerCoronaryMainModuleWidget::ShowObliqueandCurvedSlice()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	RemoveAllObliqueandCurvedSlice();
	//emit updatereformats();

	{
		vtkMRMLNode* thisaddednode;

		ObliqueSliceNode = vtkSmartPointer< vtkMRMLModelNode >::New();
		ObliqueSliceDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();

		thisaddednode = logic->GetMRMLScene()->AddNode(ObliqueSliceNode);
		addedobliqueandcurvednode.push_back(thisaddednode);
		thisaddednode = logic->GetMRMLScene()->AddNode(ObliqueSliceDisplayNode);
		addedobliqueandcurvednode.push_back(thisaddednode);
		ObliqueSliceDisplayNode->SetScene(logic->GetMRMLScene());
		ObliqueSliceNode->SetScene(logic->GetMRMLScene());
		ObliqueSliceNode->SetName("ObliqueSlice");
		ObliqueSliceNode->SetAndObserveDisplayNodeID(ObliqueSliceDisplayNode->GetID());
		ObliqueSliceNode->Modified();
		ObliqueSliceDisplayNode->Modified();

		ObliqueSliceNode->SetAndObservePolyData(ObliqueSlicePolydata);

		ObliqueSliceDisplayNode->SetVisibility(1);
		ObliqueSliceDisplayNode->SetActiveScalarName("pointcolor");
		ObliqueSliceDisplayNode->SetScalarVisibility(1);
		ObliqueSliceDisplayNode->SetAutoScalarRange(0);
		ObliqueSliceDisplayNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeGrey");
		ObliqueSliceDisplayNode->SetBackfaceCulling(0);

		//vtkSmartPointer<vtkRendererCollection> rendercollection = threeDView->GetRenderWindow()->GetRenderers();
		//vtkSmartPointer<vtkPolyDataMapper> ObliqueMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//ObliqueMapper->SetInputData(ObliqueSlicePolydata);
		//vtkSmartPointer<vtkActor> ObliqueActor = vtkSmartPointer<vtkActor>::New();
		//ObliqueActor->SetMapper(ObliqueMapper);
		//ObliqueActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
		//ObliqueActor->GetProperty()->SetLineWidth(3.0f);
		//ObliqueActor->GetProperty()->SetOpacity(0.6);
		//ObliqueActor->PickableOff();
		//vtkSmartPointer<vtkRenderer> Slicer3DRender = rendercollection->GetFirstRenderer();
		//Slicer3DRender->AddActor(ObliqueActor);
		//RenderWindowInteractorthreeD->Render();

	}

	{
		vtkMRMLNode* thisaddednode;

		CurvedSliceNode = vtkSmartPointer< vtkMRMLModelNode >::New();
		CurvedSliceDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();

		thisaddednode = logic->GetMRMLScene()->AddNode(CurvedSliceNode);
		addedobliqueandcurvednode.push_back(thisaddednode);
		thisaddednode = logic->GetMRMLScene()->AddNode(CurvedSliceDisplayNode);
		addedobliqueandcurvednode.push_back(thisaddednode);
		CurvedSliceDisplayNode->SetScene(logic->GetMRMLScene());
		CurvedSliceNode->SetScene(logic->GetMRMLScene());
		CurvedSliceNode->SetName("CurvedSlice");
		CurvedSliceNode->SetAndObserveDisplayNodeID(CurvedSliceDisplayNode->GetID());
		CurvedSliceNode->Modified();
		CurvedSliceDisplayNode->Modified();

		CurvedSliceNode->SetAndObservePolyData(CurvedSlicePolydata);

		CurvedSliceDisplayNode->SetVisibility(1);
		CurvedSliceDisplayNode->SetActiveScalarName("pointcolor");
		CurvedSliceDisplayNode->SetScalarVisibility(1);
		CurvedSliceDisplayNode->SetAutoScalarRange(0);
		CurvedSliceDisplayNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeGrey");
		CurvedSliceDisplayNode->SetBackfaceCulling(0);
	}


	return true;
}


void qSlicerCoronaryMainModuleWidget::setclcoordslot(vtkIdType pointid, double x, double y, double z)
{
//	std::cout << "in slot, pointid = " << pointid << ", newcoord = " << x << ", " << y << ", " << z << std::endl;

	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

//	double coord[3] = { x, y, z };
//	logic->centerlineModel->GetPoints()->SetPoint(pointid, coord);
	
	double coorddisplay[3] = { -x, -y, z };
	logic->centerlineModel_display->GetPoints()->SetPoint(pointid, coorddisplay);
	logic->centerlineModel_display->Modified();


/*	vtkDoubleArray *LumenRadius = vtkDoubleArray::SafeDownCast(logic->centerlineModel_display->GetPointData()->GetArray("LumenRadius"));
	vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(logic->centerlineModel_display->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(logic->centerlineModel_display->GetPointData()->GetArray("Axis2"));

	clAxis1->GetTuple(pointid, axis1);
	clAxis2->GetTuple(pointid, axis2);
	double coord[3];
*/

	return;
}

void qSlicerCoronaryMainModuleWidget::setlumenradiusslot(vtkIdType pointid, vtkIdType lumenpointid, double newradius)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	vtkDoubleArray *clArray = vtkDoubleArray::SafeDownCast(logic->centerlineModel_display->GetPointData()->GetArray("LumenRadius"));
	clArray->SetComponent(pointid, lumenpointid, newradius);
	
	vtkDoubleArray *clAxis1 = vtkDoubleArray::SafeDownCast(logic->centerlineModel->GetPointData()->GetArray("Axis1"));
	vtkDoubleArray *clAxis2 = vtkDoubleArray::SafeDownCast(logic->centerlineModel->GetPointData()->GetArray("Axis2"));
	double axis1[3], axis2[3];
	clAxis1->GetTuple(pointid, axis1);
	clAxis2->GetTuple(pointid, axis2);

	double center[3], coord[3];
	logic->centerlineModel->GetPoints()->GetPoint(pointid, center);

	double cirstep = 2 * M_PI / clArray->GetNumberOfComponents();
	for (int l = 0; l < 3; l++)
	{
		coord[l] = center[l] + newradius * (cos(lumenpointid*cirstep) * axis1[l] + sin(lumenpointid*cirstep) * axis2[l]);
	}
	
	return;
}

void qSlicerCoronaryMainModuleWidget::removemouseobserverslot()
{
//	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
//	QVTKWidget* threeDView = layoutManager->threeDWidget(0)->threeDView()->VTKWidget();
//	vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractorthreeD = threeDView->GetInteractor();

	for (int i = 0; i < this->addedvesselpickobservertag.size(); i++)
		RenderWindowInteractorthreeD->RemoveObserver(this->addedvesselpickobservertag.at(i));
	this->addedvesselpickobservertag.clear();
}


void qSlicerCoronaryMainModuleWidget::detectlumenslot(vtkIdType sid)
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (sid < 0 || sid >= logic->centerlineModel->GetNumberOfCells())
		return;

	if (logic != NULL)
	{
		if (logic->DetectLumenLogic(sid))
		{
			send_simplerendersignal();
			logic->BuildCenterlinesMeshLogic();
			SetupKeyMouseObserver();
			ShowSelectedVesselThreeD(sid);
		}
	}

}

void qSlicerCoronaryMainModuleWidget::vesseleditingwidgetclosedslot()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	this->RemoveAllSelectedVesselThreeD();
	this->RemoveAllObliqueandCurvedSlice();
	logic->BuildCenterlinesMeshLogic();
	SetupKeyMouseObserver();
}

void qSlicerCoronaryMainModuleWidget::updateobliquesliceslot(vtkPolyData* p)
{
	this->ObliqueSlicePolydata = p;
	//SavePolyData(this->ObliqueSlicePolydata, "C:\\work\\Coronary_Slicer\\testdata\\ObliqueSlicePolydata.vtp");

	//for (int i = 0; i < this->ObliqueSlicePolydata->GetPoints()->GetNumberOfPoints(); i ++)
	//{
	//	double coord[3];
	//	this->ObliqueSlicePolydata->GetPoints()->GetPoint(i, coord);
	//	coord[0] = -coord[0];
	//	coord[1] = -coord[1];
	//	this->ObliqueSlicePolydata->GetPoints()->SetPoint(i, coord);
	//}

/*	vtkMRMLNode* thisaddednode;
	 
	ObliqueSliceNode = vtkSmartPointer< vtkMRMLModelNode >::New();
	ObliqueSliceDisplayNode = vtkSmartPointer< vtkMRMLModelDisplayNode >::New();
	 
	thisaddednode =  logic->GetMRMLScene()->AddNode(ObliqueSliceNode);
	addedselectedclnode.push_back(thisaddednode);
	thisaddednode = logic->GetMRMLScene()->AddNode(ObliqueSliceDisplayNode);
	addedselectedclnode.push_back(thisaddednode);
	ObliqueSliceDisplayNode->SetScene(logic->GetMRMLScene());
	ObliqueSliceNode->SetScene(logic->GetMRMLScene());
	ObliqueSliceNode->SetName("ObliqueSlice");
	ObliqueSliceNode->SetAndObserveDisplayNodeID(ObliqueSliceDisplayNode->GetID());
	ObliqueSliceNode->Modified();
	ObliqueSliceDisplayNode->Modified();
	 
	ObliqueSliceDisplayNode->SetColor(1, 1, 0);
	ObliqueSliceDisplayNode->SetPointSize(5);
	ObliqueSliceDisplayNode->SetLineWidth(5);
	ObliqueSliceDisplayNode->SetVisibility(1);
	ObliqueSliceDisplayNode->LightingOn();
//	ObliqueSliceDisplayNode->SetRepresentation(vtkMRMLModelDisplayNode::WireframeRepresentation);
	 
	ObliqueSliceNode->SetAndObservePolyData(ObliqueSlicePolydata);
*/
	this->ObliqueSlicePolydata->Modified();

}

void qSlicerCoronaryMainModuleWidget::updatecurvedsliceslot(vtkPolyData* p)
{
	this->CurvedSlicePolydata = p;

	//for (int i = 0; i < this->CurvedSlicePolydata->GetPoints()->GetNumberOfPoints(); i++)
	//{
	//	double coord[3];
	//	this->CurvedSlicePolydata->GetPoints()->GetPoint(i, coord);
	//	coord[0] = -coord[0];
	//	coord[1] = -coord[1];
	//	this->CurvedSlicePolydata->GetPoints()->SetPoint(i, coord);
	//}

	this->CurvedSlicePolydata->Modified();
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
	//	std::cout << "Mouse click!" << std::endl;
		if (Slicer3DRenderWindowInteractor == NULL)
			return;
		if (clmodel->GetNumberOfCells() == 0)
			return;

		int pickpixel[2];
		Slicer3DRenderWindowInteractor->GetEventPosition(pickpixel);
	//	std::cout << "pickpixel = " << pickpixel[0] << ", " << pickpixel[1] << std::endl;

		vtkSmartPointer< vtkCellPicker > picker = vtkCellPicker::SafeDownCast(Slicer3DRenderWindowInteractor->GetPicker());
		picker->Pick(pickpixel[0], pickpixel[1], 0, Slicer3DRender);

		vtkIdType pickid = picker->GetCellId();
		if (pickid == -1)
			return;

		vtkSmartPointer<vtkPolyData> pickedpoly = vtkPolyData::SafeDownCast(picker->GetDataSet());
		if (pickedpoly->GetNumberOfCells() != lumenmodel->GetNumberOfCells()
			&& pickedpoly->GetNumberOfCells() != clmodel->GetNumberOfCells())
		{
			return;
		}
		if (pickedpoly->GetNumberOfCells() == lumenmodel->GetNumberOfCells())
		{
			vtkSmartPointer<vtkIdTypeArray> segmentidarray = vtkIdTypeArray::SafeDownCast(lumenmodel->GetCellData()->GetArray("SegmentId"));
			pickid = segmentidarray->GetValue(pickid);
		}
		

	//	if (pickid == LastSelectedID)
	//		return;
		LastSelectedID = pickid;

		mainwidget->send_visibilitychanged(true);

		mainwidget->SelectedVesselID = pickid;
		mainwidget->send_selectidchanged(pickid);

		mainwidget->send_clmodelchanged(clmodel);
		mainwidget->send_imagedatachanged(imagedata);
		//mainwidget->send_resetsignal();
		mainwidget->send_forcerendersignal();

		mainwidget->ShowSelectedVesselThreeD(pickid);
		mainwidget->ShowObliqueandCurvedSlice();
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


// Key Event
class vtkKeyPressedInteractionCallback : public vtkCommand
{
public:
	static vtkKeyPressedInteractionCallback *New()
	{
		return new vtkKeyPressedInteractionCallback;
	}

	vtkKeyPressedInteractionCallback()
	{
		Iren = NULL;
		VesselPicker = NULL;
	}
	~vtkKeyPressedInteractionCallback()	{}

	void SetInteractor(vtkRenderWindowInteractor *iren)
	{
		this->Iren = iren;
	}

	virtual void Execute(vtkObject *caller, unsigned long ev, void *)
	{
		if (logic->centerlineModel->GetNumberOfCells() == 0)
			return;

		if (vtkStdString(Iren->GetKeySym()) == "Control_L")
		{
			//std::cout << "Control_L is pressed!" << std::endl;
			Iren->SetPicker(VesselPicker);
			addedvesselpickobservertag->push_back(Iren->AddObserver(vtkCommand::LeftButtonPressEvent, VesselPickCallBack, 10.0f));
		}
		else if (vtkStdString(Iren->GetKeySym()) == "g" || vtkStdString(Iren->GetKeySym()) == "G")
		{
			//std::cout << "g is pressed!" << std::endl;
			mainwidget->detectlumenslot(mainwidget->SelectedVesselID);
//			mainwidget->send_visibilitychanged(true);
//			mainwidget->updatecurvedsliceslot(logic->centerlineModel);
//			mainwidget->updateobliquesliceslot();
		}
		else if (vtkStdString(Iren->GetKeySym()) == "d" || vtkStdString(Iren->GetKeySym()) == "Delete")
		{
			//std::cout << "d is pressed!" << std::endl;
			//std::cout << "SelectedVesselID = " << mainwidget->SelectedVesselID << std::endl;
			
			if (mainwidget->SelectedVesselID < 0)
				return;

			if (logic->centerlineModel && logic->centerlineModel->GetNumberOfCells() > 0)
			{
				if (mainwidget->SelectedVesselID < logic->centerlineModel->GetNumberOfCells())
				{
					logic->centerlineModel->BuildCells();
					logic->DeleteCenterlineOneSegmentLogic(mainwidget->SelectedVesselID);

			//		SavePolyData(logic->centerlineModel, "C:\\work\\Coronary_Slicer\\testdata\\logic_centerlineModel1.vtp");

					if (logic->centerlineModel->GetCellData()->HasArray("SegmentId"))
						logic->centerlineModel->GetCellData()->RemoveArray("SegmentId");

					logic->centerlineId = vtkSmartPointer<vtkIdFilter>::New();
					logic->centerlineId->SetInputData(logic->centerlineModel);
					logic->centerlineId->PointIdsOff();
					logic->centerlineId->CellIdsOn();
					logic->centerlineId->FieldDataOn();
					logic->centerlineId->SetIdsArrayName("SegmentId");
					logic->centerlineId->Update();
					logic->centerlineModel = vtkPolyData::SafeDownCast(logic->centerlineId->GetOutput());
					
					logic->AddCircumParamtoClModel();
					logic->AddLongiParamtoClModel();
					
					logic->BuildCenterlinesMeshLogic();
					mainwidget->SetupKeyMouseObserver();

					mainwidget->SelectedVesselID = -1;
					mainwidget->send_visibilitychanged(false);

				}
			}
		}
	}

private:
	vtkRenderWindowInteractor* Iren;
public:
	vtkCellPicker* VesselPicker;
	CVesselPickCallBack* VesselPickCallBack;
	vector<unsigned long>* addedvesselpickobservertag;

	qSlicerCoronaryMainModuleWidget* mainwidget;
	vtkSlicerCoronaryMainLogic *logic;
	
//	vtkPolyData* clmodel;
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
			//std::cout << "Control_L is released!" << std::endl;
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
::RemoveKeyMouseObserver()
{
//	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
//	QVTKWidget* threeDView = layoutManager->threeDWidget(0)->threeDView()->VTKWidget();
//	vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractorthreeD = threeDView->GetInteractor();

	for (int i = 0; i < this->addedvesselpickobservertag.size(); i++)
		RenderWindowInteractorthreeD->RemoveObserver(this->addedvesselpickobservertag.at(i));
	this->addedvesselpickobservertag.clear();

	for (int i = 0; i < this->addedctrlobservertag.size(); i++)
		RenderWindowInteractorthreeD->RemoveObserver(this->addedctrlobservertag.at(i));
	this->addedctrlobservertag.clear();
	return true;
}


bool qSlicerCoronaryMainModuleWidget
::SetupKeyMouseObserver()
{
	Q_D(qSlicerCoronaryMainModuleWidget);
	vtkSlicerCoronaryMainLogic *logic = d->logic();

	if (logic->centerlineTube->GetOutput(0)->GetNumberOfCells() == 0)
		return false;

	qSlicerCoronaryMainModuleWidget::RemoveAllSelectedVesselThreeD();
//	qSlicerCoronaryMainModuleWidget::RemoveAllObliqueandCurvedSlice();
	qSlicerCoronaryMainModuleWidget::RemoveKeyMouseObserver();

	// set ctrl observer
//	qSlicerLayoutManager* layoutManager = qSlicerApplication::application()->layoutManager();
//	QVTKWidget* threeDView = layoutManager->threeDWidget(0)->threeDView()->VTKWidget();
//	vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractorthreeD = threeDView->GetInteractor();
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
	
	

	VesselPicker = vtkSmartPointer<vtkCellPicker>::New();
	VesselPicker->SetTolerance(0.005);
	VesselPicker->PickClippingPlanesOff();

	VesselPickCallBack = vtkSmartPointer<CVesselPickCallBack>::New();
	VesselPickCallBack->Slicer3DRenderWindowInteractor = RenderWindowInteractorthreeD;
	VesselPickCallBack->Slicer3DRender = rendercollection->GetFirstRenderer();
	VesselPickCallBack->mainwidget = this;
	VesselPickCallBack->clmodel = logic->centerlineModel;
	VesselPickCallBack->lumenmodel = logic->centerlineTube->GetOutput(2);
	VesselPickCallBack->imagedata = logic->imageData;


	vtkSmartPointer<vtkKeyPressedInteractionCallback> KeyPressedInteractionCallback = vtkSmartPointer<vtkKeyPressedInteractionCallback>::New();
	KeyPressedInteractionCallback->SetInteractor(RenderWindowInteractorthreeD);
	KeyPressedInteractionCallback->mainwidget = this;
	KeyPressedInteractionCallback->logic = logic;
	KeyPressedInteractionCallback->VesselPicker = VesselPicker;
	KeyPressedInteractionCallback->VesselPickCallBack = VesselPickCallBack;
	KeyPressedInteractionCallback->addedvesselpickobservertag = &addedvesselpickobservertag;
//	KeyPressedInteractionCallback->clmodel = logic->centerlineModel;
	addedctrlobservertag.push_back(RenderWindowInteractorthreeD->AddObserver(vtkCommand::KeyPressEvent, KeyPressedInteractionCallback));

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

	
//	d->logic()->BuildCenterlinesMeshLogic();
//	SetupKeyMouseObserver();

//	vtkDoubleArray* temp = vtkDoubleArray::SafeDownCast(d->logic()->centerlineModel->GetPointData()->GetArray("LongiParam"));
//	std::cout << "temp.num = " << temp->GetNumberOfValues() << std::endl;
//	std::cout << "point.num = " << d->logic()->centerlineModel->GetPoints()->GetNumberOfPoints() << std::endl;
	

//	SavePolyData(d->logic()->centerlineTube->GetOutput(0), "C:\\work\\Coronary_Slicer\\testdata\\centerlineTube0.vtp");
//	d->logic()->TestLogic();

	std::cout << "TestButtonFunc end" << std::endl;

	return true;
}
