#ifndef GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iterator_range.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class InputRange, class ElementMap>
            class Generalized_region_growing_point_traits_2 {
            public:
                typedef Kernel::Point_2 Element;
                typedef InputRange Input_range;
                typedef ElementMap Element_map;
            }
        }
    }
}

#endif