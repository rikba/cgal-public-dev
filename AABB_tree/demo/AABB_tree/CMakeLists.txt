# This is the CMake script for compiling the AABB tree demo.

project( AABB_demo )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )
# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED COMPONENTS Qt5 )
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
find_package(Qt5 COMPONENTS Xml Script OpenGL Gui)

# Find OpenGL
find_package(OpenGL)

macro(ADD_QT_PROJECT target_name source_files qt_resources resource_files lib_files include_directories out_file_path)
    #ADDING CONFIG AND TARGET NAME
    file(WRITE ${out_file_path} "CONFIG += qt c++11 resources\n" )
    file(APPEND ${out_file_path} "TARGET = ${target_name}\n")
    file(APPEND ${out_file_path} "QT += ${qt_resources}\n")


    #SETTING UP TARGET TYPE
    get_target_property(TARGET_CMAKE_TYPE ${target_name} TYPE)
    if(${TARGET_CMAKE_TYPE} STREQUAL "EXECUTABLE")
    file(APPEND ${out_file_path} "TEMPLATE = app\n")
    else()
    file(APPEND ${out_file_path} "TEMPLATE = lib\n")
    endif()

    #ADDING SOURCES
    file(APPEND ${out_file_path} "SOURCES += ${source_files}\n")

    #ADDING RESOURCES
    file(APPEND ${out_file_path} "RESOURCES += ${resource_files}\n")

    #ADDING LIBRARIES
    file(APPEND ${out_file_path} "LIBS += " "${lib_files}" "\n")
    #ADDING INCLUDE DIRECTORIES
    file(APPEND ${out_file_path} "INCLUDEPATH =  ${include_directories} \n")

endmacro(ADD_QT_PROJECT)

# Find QGLViewer
if(Qt5_FOUND)
  find_package(QGLViewer )
endif(Qt5_FOUND)

if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
  set (QT_RESOURCE "qml quick widgets gui opengl xml script")
  set (SOURCES "AABB_demo.cpp ")
  set(RESOURCE_FILES2 "AABB_demo.qrc" )
  set( OUT_PATH "AABB_demo.pro")
  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt5_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt5_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt5_add_resources ( RESOURCE_FILES AABB_demo.qrc )

  add_file_dependencies( AABB_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )

  add_executable  ( AABB_demo AABB_demo.cpp ${UI_FILES} ${RESOURCE_FILES} )
  get_property(LIB_FILES TARGET AABB_demo PROPERTY LINK_LIBRARIES)
  get_property(INCLUDE_DIRECTORY TARGET AABB_demo PROPERTY INCLUDE_DIRECTORIES )
  ADD_QT_PROJECT( AABB_demo ${SOURCES} ${QT_RESOURCE} ${RESOURCE_FILES2} "${LIB_FILES}" "${INCLUDE_DIRECTORY}" ${OUT_PATH})
  qt5_use_modules(AABB_demo Gui OpenGL Xml Script )
  # Link with Qt libraries
  target_link_libraries( AABB_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( AABB_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( AABB_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS AABB_demo )


else (CGAL_FOUND AND CGAL_CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(AABB_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(AABB_MISSING_DEPS "the CGAL library, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt5_FOUND)
    set(AABB_MISSING_DEPS "the CGAL Qt5 library, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(AABB_MISSING_DEPS "Qt5, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(AABB_MISSING_DEPS "OpenGL, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(AABB_MISSING_DEPS "QGLViewer, ${AABB_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${AABB_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
