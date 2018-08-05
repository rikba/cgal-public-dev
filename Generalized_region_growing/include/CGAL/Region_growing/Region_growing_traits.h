#ifndef CGAL_REGION_GROWING_TRAITS_H
#define CGAL_REGION_GROWING_TRAITS_H

namespace CGAL {
    namespace Region_growing {

            /*!
                \ingroup PkgGeneralizedRegionGrowing
                \brief Traits class describing the input values and the primary type of geometric element on which we grow regions.
                \tparam InputRange
                \tparam ElementMap
                \tparam Kernel_
            */

            template<class InputRange, class ElementMap, class Kernel_>
            class Region_growing_traits {
            public:
                
                using Kernel                  = Kernel_;
                ///< The kernel on which the algorithms working.

                using Input_range             = InputRange;
                ///< An `Iterator_range` of bidirectional iterator.

                using Element_map             = ElementMap;
                ///< A model of `LvaluePropertyMap` that maps to a geometric type.

                using Element                 = typename Element_map::value_type;
                ///< The value type of the iterators in `Input_range`.
                
                // ElementMap::value_type must be K::Point_2 or K::Point_3 or Face_index
                // ElementMap::key_type must be the same as std::iterator_traits<InputRange::iterator>::value_type
            };

    }

}

#endif