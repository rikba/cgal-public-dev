# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Alpha_shapes_2_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 COMPONENTS Xml Script OpenGL)
set (QT_RESOURCE "qml quick widgets script xml opengl")
set (SOURCES "Alpha_shapes_2.cpp Alpha_shapes_2.moc")
set(RESOURCE_FILES "Alpha_shapes_2.qrc" )
set( OUT_PATH "Alpha_shapes_2.pro")
#set(LIB_FILES ${QT_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

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

include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )


#--------------------------------
# The "Delaunay" demo: Alpha_shapes_2
#--------------------------------
# UI files (Qt Designer files)
qt5_wrap_ui( DT_UI_FILES Alpha_shapes_2.ui )

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( DT_RESOURCE_FILES ./Alpha_shapes_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt5_generate_moc( Alpha_shapes_2.cpp Alpha_shapes_2.moc )

# The executable itself.
add_executable  ( Alpha_shapes_2 Alpha_shapes_2.cpp Alpha_shapes_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

qt5_use_modules(Alpha_shapes_2 Xml Script OpenGL)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Alpha_shapes_2 )

# Link with Qt libraries
target_link_libraries( Alpha_shapes_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Alpha_shapes_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
get_property(LIB_FILES TARGET Alpha_shapes_2 PROPERTY LINK_LIBRARIES)
get_property(INCLUDE_DIRECTORY TARGET Alpha_shapes_2 PROPERTY INCLUDE_DIRECTORIES )
ADD_QT_PROJECT( Alpha_shapes_2 ${SOURCES} ${QT_RESOURCE} ${RESOURCE_FILES} "${LIB_FILES}" "${INCLUDE_DIRECTORY}" ${OUT_PATH})
else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
