#pragma once
#pragma warning(disable:4005)

/*
#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif
*/

namespace JPTD {

typedef unsigned char uchar;
typedef unsigned int uint;

typedef enum {
	PLUS = 1,
	ZERO = 0,
	MINUS = -1
} Sign;

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

/*#ifndef COS_PI_360
#define COS_PI_360 0.99996192306
#endif*/

/*
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point;
typedef bg::model::box<Boost_Point> Boost_Box;
typedef std::pair<Boost_Box, uint> Boost_Value;
typedef bgi::rtree<Boost_Value, bgi::quadratic<16> > Boost_RTree;
*/

#ifdef KINETIC_PARTITION_EXPORTS
#define KINETIC_PARTITION_API __declspec(dllexport)
#else
#define KINETIC_PARTITION_API __declspec(dllimport)
#endif

}