#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/IO/Color.h>

namespace JPTD {

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Lazy_exact_nt<CGAL::Epeck_ft> FT;

typedef K::Point_2 CGAL_Point_2;
typedef K::Point_3 CGAL_Point_3;
typedef K::Segment_2 CGAL_Segment_2;
typedef K::Line_2 CGAL_Line_2;
typedef K::Line_3 CGAL_Line_3;
typedef K::Vector_2 CGAL_Vector_2;
typedef K::Vector_3 CGAL_Vector_3;
typedef K::Plane_3 CGAL_Plane;
typedef K::Direction_2 CGAL_Direction_2;

typedef CGAL::Color CGAL_Color;

typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 CGAL_Inexact_Point_3;

typedef CGAL::Polygon_2<K> CGAL_Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> CGAL_Polygon_with_holes_2;

}