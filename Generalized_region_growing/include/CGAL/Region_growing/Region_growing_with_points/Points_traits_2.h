#ifndef GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class InputRange, class ElementMap, class Kernel>
            class Generalized_region_growing_point_traits_2 {
            public:
                using Kernel      = Kernel;
                using Element     = Kernel::Point_2;
                using Input_range = InputRange;
                // ElementMap::value_type must be Point_2 (Element)
                using Element_map = ElementMap;
            }
        }
    }
}

#endif