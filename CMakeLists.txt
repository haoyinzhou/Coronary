cmake_minimum_required(VERSION 2.8.9)

project(Coronary)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/Coronary")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Haoyin Zhou (BWH)")
set(EXTENSION_DESCRIPTION "CTA imaging -based Coronary segmentation and visulazation ")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/Coronary.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/Coronary/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

find_package( OpenCV REQUIRED )


#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CoronaryMain)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
