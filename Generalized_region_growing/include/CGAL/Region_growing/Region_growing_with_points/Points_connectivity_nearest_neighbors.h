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

            template<class Traits>
            class Points_connectivity_nearest_neighbors {

            public:
                using Input_range             = typename Traits::Input_range;
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;
                using Element                 = typename Traits::Element;

                using Element_with_properties = typename Element_map::key_type;

                using Neighbors               = std::vector<Element_with_properties>;
                using Neighbor_range          = CGAL::Iterator_range<typename Neighbors::const_iterator>;

                using Search_base             = typename std::conditional<std::is_same<typename Kernel::Point_2, Element>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;

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

                // Element == Point_3 or Point_2
                using Point             = typename Traits::Element;
                using FT                = typename Kernel::FT;
                using Kd_tree_config    = Search_structures<Point, Element_with_properties>;
                using Neighbor_search   = typename Kd_tree_config::Neighbor_search;
                using Tree              = typename Kd_tree_config::Tree;

                Points_connectivity_nearest_neighbors(const Input_range &input_range, const unsigned int number_of_neighbors) :
                        m_input_range(input_range),
                        m_number_of_neighbors(number_of_neighbors) {
                    m_tree.insert(m_input_range.begin(), m_input_range.end());
                }

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
