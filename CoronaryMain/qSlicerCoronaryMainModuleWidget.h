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

class qSlicerCoronaryMainModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_CORONARYMAIN_EXPORT qSlicerCoronaryMainModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerCoronaryMainModuleWidget(QWidget *parent=0);
  virtual ~qSlicerCoronaryMainModuleWidget();

public slots:
	void TestButtonFuc();

protected:
  QScopedPointer<qSlicerCoronaryMainModuleWidgetPrivate> d_ptr;
  virtual void setup();



private:
  Q_DECLARE_PRIVATE(qSlicerCoronaryMainModuleWidget);
  Q_DISABLE_COPY(qSlicerCoronaryMainModuleWidget);
};

#endif
