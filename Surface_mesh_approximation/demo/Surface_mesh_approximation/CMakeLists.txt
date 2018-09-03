# This is the CMake script for compiling the VSA demo.

project( Surface_mesh_approximation_Demo )

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

add_definitions(-std=c++11)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
find_package(Qt5 REQUIRED QUIET COMPONENTS Xml OpenGL)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
if(Qt5_FOUND)
  find_package(QGLViewer REQUIRED)
endif(Qt5_FOUND)

# Active verbose output ? (turned ON by default)
option(ACTIVATE_VERBOSE_OUTPUT
  "Enable verbose output"
  ON)

if(ACTIVATE_VERBOSE_OUTPUT)
  add_definitions(-DCGAL_SURFACE_MESH_APPROXIMATION_DEBUG)
endif()

# Activate concurrency ? (turned OFF by default)
option(ACTIVATE_CONCURRENCY
  "Enable concurrency"
  OFF)

if(ACTIVATE_CONCURRENCY)
  add_definitions(-DCGAL_LINKED_WITH_TBB)
  find_package(TBB REQUIRED)
endif()

if(TBB_FOUND)
  include(${TBB_USE_FILE})
endif()

if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  # Include this package's headers first
  include_directories( BEFORE ../../include ./ )

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  # ui files, created with Qt Designer
  qt5_wrap_ui( UI_FILES MainWindow.ui dialSettings.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES VSA_demo.qrc )

  # cpp files
  add_executable( VSA_demo
    VSA_demo.cpp
    MainWindow.cpp Viewer.cpp
    dialSettings.cpp
    Scene.cpp
    ColorCheatSheet.cpp
    ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

  qt5_use_modules( VSA_demo Xml OpenGL )

  # Link with Qt libraries
  target_link_libraries( VSA_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( VSA_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( VSA_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

  # Link with TBB
  if(ACTIVATE_CONCURRENCY)
    if(TBB_FOUND)
      target_link_libraries(VSA_demo ${TBB_LIBRARIES})
    endif()
  endif()

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS VSA_demo )

else (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(VSA_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(VSA_MISSING_DEPS "the CGAL library, ${VSA_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt5_FOUND)
    set(VSA_MISSING_DEPS "CGAL_Qt5, ${VSA_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(VSA_MISSING_DEPS "Qt5, ${VSA_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(VSA_MISSING_DEPS "OpenGL, ${VSA_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(VSA_MISSING_DEPS "QGLViewer, ${VSA_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${VSA_MISSING_DEPS} and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
