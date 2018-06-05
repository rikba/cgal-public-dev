#ifndef GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class InputRange, class ElementMap, class Kernel>
            class Generalized_region_growing_point_traits_2 {
            public:
                using Kernel                  = Kernel;
                using Input_range             = InputRange;
                using Element_map             = ElementMap;
                using Element_with_properties = Element_map::key_type;
                using Element                 = Kernel::Point_2;
                // ElementMap::value_type must be the same as Element
                // ElementMap::key_type must be the same as std::iterator_traits<InputRange::iterator>::value_type
            }
        }
    }
}

#endif