#pragma once
#pragma warning(disable:4005)

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/IO/Color.h>
#include <eigen3/Eigen/Dense>

#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

namespace JPTD {

typedef unsigned char uchar;
typedef unsigned int uint;

#ifndef jmin
#define jmin(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef jmax
#define jmax(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef jclamp
#define jclamp(a, x, b) ((x) < (a) ? (a) : ((x) > (b) ? (b) : (x)))
#endif
#ifndef jin
#define jin(a, x, b) ((a) <= (x) && (x) <= (b))
#endif
#ifndef PI
#define PI 3.141592653589783238462643383279
#endif

#ifndef COS_PI_360
#define COS_PI_360 0.99996192306
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Lazy_exact_nt<CGAL::Epeck_ft> FT;

typedef K::Point_2 CGAL_Point_2;
typedef K::Point_3 CGAL_Point_3;
typedef K::Segment_2 CGAL_Segment_2;
typedef K::Segment_3 CGAL_Segment_3;
typedef K::Line_2 CGAL_Line_2;
typedef K::Line_3 CGAL_Line_3;
typedef K::Vector_2 CGAL_Vector_2;
typedef K::Vector_3 CGAL_Vector_3;
typedef K::Plane_3 CGAL_Plane;
typedef K::Triangle_2 CGAL_Triangle_2;
typedef K::Triangle_3 CGAL_Triangle_3;
typedef CGAL::Color CGAL_Color;

typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 CGAL_EPICK_Point_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 CGAL_EPICK_Point_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel::Vector_2 CGAL_EPICK_Vector_2;

typedef CGAL::Polygon_2<K> CGAL_Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> CGAL_Polygon_with_holes_2;

typedef Eigen::Matrix3d Matrix_3;

typedef enum {
	PLUS = 1,
	ZERO = 0,
	MINUS = -1
} Sign;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point;
typedef bg::model::box<Boost_Point> Boost_Box;
typedef std::pair<Boost_Box, uint> Boost_Value;
typedef bgi::rtree<Boost_Value, bgi::quadratic<16> > Boost_RTree;


#ifdef KINETIC_PARTITION_EXPORTS
#define KINETIC_PARTITION_API __declspec(dllexport)
#else
#define KINETIC_PARTITION_API __declspec(dllimport)
#endif

}