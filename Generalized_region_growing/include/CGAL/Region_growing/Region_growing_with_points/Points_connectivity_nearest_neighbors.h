#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_NEAREST_NEIGHBORS_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_NEAREST_NEIGHBORS_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            /*!
                \ingroup PkgGeneralizedRegionGrowingPoints
                \brief K nearest neighbors (kNN) searching method on a set of `Point_2` or `Point_3`
                \tparam Traits_ `CGAL::Region_growing::Region_growing_traits`
                \cgalModels `RegionGrowingConnectivity`
            */

            template<class Traits_>
            class Points_connectivity_nearest_neighbors {

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
                        using Neighbor_search       = CGAL::Orthogonal_k_neighbor_search<Search_traits_adapter>;
                        using Tree                  = typename Neighbor_search::Tree;
                    };

                    // Partial specialization when PointType and ElementWithProperties are the same
                    template<class PointType>
                    struct Search_structures<PointType, PointType> {
                        using Neighbor_search       = CGAL::Orthogonal_k_neighbor_search<Search_base>;
                        using Tree                  = typename Neighbor_search::Tree;
                    };

                #endif

                using FT                = typename Kernel::FT; ///< Number type

                using Kd_tree_config    = Search_structures<Point, Element_with_properties>;
                ///< Kd tree configuration

                using Neighbor_search   = typename Kd_tree_config::Neighbor_search;
                ///< A search class CGAL::Orthogonal_k_neighbor_search that implements a kd tree

                using Tree              = typename Kd_tree_config::Tree;
                ///< The kd tree member type of the search class

                /*!
                    The constructor initializes a kd tree using user input and stores `number_of_neighbors` (the value of "k", as in "kNN") for later use.
                */
                Points_connectivity_nearest_neighbors(const Input_range &input_range, const unsigned int number_of_neighbors) :
                m_input_range(input_range),
                m_number_of_neighbors(number_of_neighbors) {
                    m_tree.insert(m_input_range.begin(), m_input_range.end());
                }

                /*!
                    The function takes a query point and return k closest points around that. The result is stored in `neighbors`.
                    \tparam Neighbors_ CGAL::Region_growing::Generalized_region_growing::Neighbors
                */
                template < class Neighbors_ >
                void get_neighbors(const Element_with_properties &center, Neighbors_& neighbors) {

                    neighbors.clear();
                    Neighbor_search search(m_tree, get(m_elem_map, center), m_number_of_neighbors);
                    for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
                        neighbors.push_back(it->first);

                }

            private:
                const Element_map               m_elem_map = Element_map();
                const Input_range &             m_input_range;
                Tree                            m_tree;
                const unsigned int              m_number_of_neighbors;
            };

        } // namespace Region_growing_with_points 

    } // namespace Region_growing

} // namespace CGAL

#endif
