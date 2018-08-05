#ifndef CGAL_GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H
#define CGAL_GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            /*!
                \ingroup PkgGeneralizedRegionGrowingPoints
                \brief Circular neighbor searching method on a set of `Point_2` or `Point_3`
                \tparam Traits_ `CGAL::Region_growing::Region_growing_traits`
                \cgalModels `RegionGrowingConnectivity`
            */

            template<class Traits_>
            class Points_connectivity_circular_query {

            public:
                using Input_range             = typename Traits_::Input_range;
                using Kernel                  = typename Traits_::Kernel;
                using Element_map             = typename Traits_::Element_map;

                using Point                   = typename Traits_::Element;
                ///< Point type, can only be `Point_2` or `Point_3`

                using Element_with_properties = typename Element_map::key_type;
                ///< The value type of iterators in `Input_range`

                #ifdef DOXYGEN_RUNNING
                    using Search_base         = unspecified_type;
                    ///< Can be `CGAL::Search_traits_2` or `CGAL::Search_traits_3`, automatically deduced based on whether the point type is Point_3 or Point_2.

                    using Search_structures   = unspecified_type;
                    ///< Kd tree configuration class, automatically deduced based on whether `Element_map` is a `CGAL::Identity_property_map` or not. If the `Element_map` is an identity map, it will use the traits class Search_base directly to construct the tree, otherwise the program will create a traits adapter and use that instead.
                #else
                    using Search_base             = typename std::conditional<std::is_same<typename Kernel::Point_2, Point>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;
                    // Primary template
                    template<class PointType, class ElementWithProperties>
                    struct Search_structures {
                        using Search_traits_adapter = CGAL::Search_traits_adapter<Element_with_properties, Element_map, Search_base>;
                        using Fuzzy_circle          = CGAL::Fuzzy_sphere<Search_traits_adapter>;
                        using Tree                  = CGAL::Kd_tree<Search_traits_adapter>;
                    };

                    // Partial specialization when PointType and ElementWithProperties are the same
                    template<class PointType>
                    struct Search_structures<PointType, PointType> {
                        using Fuzzy_circle  = CGAL::Fuzzy_sphere<Search_base>;
                        using Tree          = CGAL::Kd_tree<Search_base>;
                    };
                #endif

                // Element == Point_3 or Point_2

                using FT                = typename Kernel::FT; ///< Number type

                using Kd_tree_config    = Search_structures<Point, Element_with_properties>;
                ///< Kd tree configuration
                
                using Fuzzy_circle      = typename Kd_tree_config::Fuzzy_circle;
                ///< A `CGAL::Fuzzy_sphere` representing the search space of the algorithm.
                
                using Tree              = typename Kd_tree_config::Tree;
                ///< A `CGAL::Kd_tree` holding the points given in the input range.

                /*!
                    The constructor takes the point set given in `input_range` and the searching radius, then initializes a kd tree upon the point set.
                */
                Points_connectivity_circular_query(const Input_range &input_range, const FT radius) :
                m_input_range(input_range),
                m_radius(radius) {
                    m_tree.insert(m_input_range.begin(), m_input_range.end());
                }

                /*!
                    From a query point `center`, this function creates a `CGAL::Fuzzy_sphere` using the radius previously given in the constructor. It then uses CGAL::Kd_tree::search() to look for the neighbors and push them to `neighbors`.
                    \tparam Neighbors_ CGAL::Region_growing::Generalized_region_growing::Neighbors
                */
                template < class Neighbors_ >
                void get_neighbors(const Element_with_properties &center, Neighbors_& neighbors) {

                    neighbors.clear();
                    Fuzzy_circle circle(center, m_radius);
                    m_tree.search(std::back_inserter(neighbors), circle);
                    
                }

            private:
                const Element_map               m_elem_map = Element_map();
                const Input_range &             m_input_range;
                Tree                            m_tree;
                const FT                        m_radius;
            };

        } // namespace Region_growing_with_points 

    } // namespace Region_growing

} // namespace CGAL

#endif