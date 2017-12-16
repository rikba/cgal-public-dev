# This is the CMake script for compiling the VSA demo.

project( Surface_mesh_approximation_Demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Activate concurrency?
option(ACTIVATE_CONCURRENCY
  "Enable concurrency"
  ON)

# Find TBB
# Tested version: tbb2018_20171205oss
# Find TBB by setting $ENV{TBB_DIR} to the directory where TBBConfig.cmake resides.
# Placed before `include(${CGAL_USE_FILE})` to avoid using the CGAL Installation/cmake/modules/FindTBB.cmake module.
# Reference: $ENV{TBB_DIR}/README.rst
# find_package will search for available TBBConfig using variables CMAKE_PREFIX_PATH and TBB_DIR.
if(ACTIVATE_CONCURRENCY)
  # add definition required by CGAL
  # https://github.com/CGAL/cgal/wiki/Concurrency-in-CGAL#which-mechanism-to-enable-or-disable-parallelism
  add_definitions(-DCGAL_LINKED_WITH_TBB)
  find_package(TBB REQUIRED tbb)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

include( ${CGAL_USE_FILE} )

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(Qt5_FOUND)
  find_package(QGLViewer)
endif(Qt5_FOUND)

if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt5_wrap_ui( UI_FILES MainWindow.ui dialSettings.ui)

  include(AddFileDependencies)

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt5_generate_moc( "dialSettings.h" "${CMAKE_CURRENT_BINARY_DIR}/dialSettings_moc.cpp" )
  add_file_dependencies( dialSettings_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/dialSettings.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES VSA_demo.qrc )

  add_file_dependencies( VSA_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/dialSettings_moc.cpp")

  add_executable  ( VSA_demo VSA_demo.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})
  qt5_use_modules(VSA_demo Xml Script OpenGL)
  # Link with Qt libraries
  target_link_libraries( VSA_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( VSA_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( VSA_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

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

# Reference: $ENV{TBB_DIR}/README.rst
# Link Intel TBB imported targets to the executable;
# "TBB::tbb" can be used instead of "${TBB_IMPORTED_TARGETS}".
if(ACTIVATE_CONCURRENCY)
  if(TBB_FOUND)
    target_link_libraries(VSA_demo ${TBB_IMPORTED_TARGETS})
  endif()
endif()
