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
#include <QtPlugin>

// CoronaryMain Logic includes
#include <vtkSlicerCoronaryMainLogic.h>

// CoronaryMain includes
#include "qSlicerCoronaryMainModule.h"
#include "qSlicerCoronaryMainModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerCoronaryMainModule, qSlicerCoronaryMainModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerCoronaryMainModulePrivate
{
public:
  qSlicerCoronaryMainModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerCoronaryMainModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModulePrivate::qSlicerCoronaryMainModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerCoronaryMainModule methods

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModule::qSlicerCoronaryMainModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerCoronaryMainModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerCoronaryMainModule::~qSlicerCoronaryMainModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerCoronaryMainModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerCoronaryMainModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerCoronaryMainModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerCoronaryMainModule::icon() const
{
  return QIcon(":/Icons/CoronaryMain.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerCoronaryMainModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerCoronaryMainModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerCoronaryMainModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerCoronaryMainModule
::createWidgetRepresentation()
{
  return new qSlicerCoronaryMainModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerCoronaryMainModule::createLogic()
{
  return vtkSlicerCoronaryMainLogic::New();
}
