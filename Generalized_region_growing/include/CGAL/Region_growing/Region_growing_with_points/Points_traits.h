#ifndef CGAL_REGION_GROWING_POINTS_TRAITS_H
#define CGAL_REGION_GROWING_POINTS_TRAITS_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {

            template<class InputRange, class ElementMap, class K>
            class Points_traits {
            public:
                using Kernel                  = K;
                using Input_range             = InputRange;
                using Element_map             = ElementMap;
                using Element                 = typename Element_map::value_type;
                
                // ElementMap::value_type must be K::Point_2 or K::Point_3
                // ElementMap::key_type must be the same as std::iterator_traits<InputRange::iterator>::value_type
            };

        }
    }
}

#endif