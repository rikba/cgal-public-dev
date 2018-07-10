# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project(Generalized_region_growing_Examples)

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)

if (CGAL_FOUND)

    include(${CGAL_USE_FILE})

    include(CGAL_CreateSingleSourceCGALProgram)

    include_directories(BEFORE "../../include")

    include_directories(BEFORE "include")

    create_single_source_cgal_program("Simple_example.cpp")
    create_single_source_cgal_program("Region_growing_with_points_2.cpp")
    create_single_source_cgal_program("Region_growing_with_points_3.cpp")
    create_single_source_cgal_program("Region_growing_with_meshes.cpp")

else ()

    message(WARNING "This program requires the CGAL library, and will not be compiled.")

endif ()

