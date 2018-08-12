project(Generalized_region_growing_Benchmark)

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)

if (CGAL_FOUND)

    include(${CGAL_USE_FILE})

    include(CGAL_CreateSingleSourceCGALProgram)

    include_directories(BEFORE "../../include")

    create_single_source_cgal_program("benchmark_region_growing_with_points_2.cpp")

else ()

    message(WARNING "This program requires the CGAL library, and will not be compiled.")

endif ()

