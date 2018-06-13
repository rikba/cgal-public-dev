#ifndef GENERALIZED_REGION_GROWING_POINTS_TRAITS_3_H
#define GENERALIZED_REGION_GROWING_POINTS_TRAITS_3_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template<class InputRange, class ElementMap, class K>
            class Points_traits_3 {
            public:
                using Kernel                  = K;
                using Input_range             = InputRange;
                using Element_map             = ElementMap;
                using Element                 = typename Kernel::Point_3;
                // ElementMap::value_type must be the same as Element
                // ElementMap::key_type must be the same as std::iterator_traits<InputRange::iterator>::value_type
            };
        }
    }
}

#endif