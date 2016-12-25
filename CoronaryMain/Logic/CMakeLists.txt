project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  Common.cxx
  Common.h
  LearningImpl.cxx
  LearningImpl.h
  ExtendSplineFilter.h
  ExtendSplineFilter.cxx
  itkBinaryThinningImageFilter3D.h
  itkBinaryThinningImageFilter3D.txx
  ExtendTubeFilter.h
  ExtendTubeFilter.cxx
  ImageStretchCurvedReformat.h
  ImageStretchCurvedReformat.cxx
  ImageCurvedReformat.h
  ImageCurvedReformat.cxx
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${OpenCV_LIBS}
  )

find_package( OpenCV REQUIRED )


#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )