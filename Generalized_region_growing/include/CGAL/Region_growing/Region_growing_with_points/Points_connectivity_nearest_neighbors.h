#ifndef CGAL_REGION_GROWING_POINTS_CONNECTIVITY_NEAREST_NEIGHBORS_H
#define CGAL_REGION_GROWING_POINTS_CONNECTIVITY_NEAREST_NEIGHBORS_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/property_map.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            template<class Traits>
            class Points_connectivity_nearest_neighbors {

            public:
                using Input_range             = typename Traits::Input_range;
                using Kernel                  = typename Traits::Kernel;
                using Element                 = typename Traits::Element; // Point_2 or Point_3
                using Element_map_original    = typename Traits::Element_map;
                using Element_index           = size_t;
                
                struct Element_map {
                    using key_type = Element_index;
                    using value_type = Element;
                    using reference = const value_type&;
                    using category = boost::lvalue_property_map_tag;
                    using Self = Element_map;

                    Input_range input_range;
                    Element_map_original elem_map_orig = Element_map_original();

                    Element_map(const Input_range& ir) : input_range(ir) { }

                    value_type& operator[](key_type& k) const { return elem_map_orig[*(input_range.begin() + k)]; }

                    friend reference get(const Self& map, const key_type& k) {
                        return get(map.elem_map_orig, *(map.input_range.begin() + k));
                    }
                };

                using Search_base             = typename std::conditional<std::is_same<typename Kernel::Point_2, Element>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;
                using Search_traits_adapter   = CGAL::Search_traits_adapter<int, Element_map, Search_base>;
                using Distance_base           = CGAL::Euclidean_distance<Search_base>;
                using Distance_adapter        = CGAL::Distance_adapter<Element_index, Element_map, Distance_base>;
                using Tree                    = CGAL::Kd_tree<Search_traits_adapter>;
                using Splitter                = CGAL::Sliding_midpoint<Search_traits_adapter>;
                using Neighbor_search         = CGAL::Orthogonal_k_neighbor_search<Search_traits_adapter, Distance_adapter, Splitter, Tree>;

                using FT                      = typename Kernel::FT;

                Points_connectivity_nearest_neighbors(const Input_range &input_range, const unsigned int number_of_neighbors) :
                m_input_range(input_range),
                m_number_of_neighbors(number_of_neighbors),
                m_elem_map(Element_map(input_range)),
                m_tree(Splitter(), Search_traits_adapter(Element_map(input_range), Search_base())),
                m_distance_adapter(Element_map(input_range)){

                    size_t i = 0;
                    for (typename Input_range::iterator it = m_input_range.begin(); it != m_input_range.end(); ++it, ++i)
                        m_indices.push_back(i);

                    m_tree.insert(m_indices.begin(), m_indices.end());

                }

                template < class Neighbors_ >
                void get_neighbors(const Element_index center, Neighbors_& neighbors) {

                    neighbors.clear();
                    Neighbor_search search(m_tree, get(m_elem_map, center), m_number_of_neighbors, 0.0, true, m_distance_adapter, true);
                    for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
                        neighbors.push_back(it->first);

                }

            private:
                const Element_map               m_elem_map;
                const Input_range &             m_input_range;
                Tree                            m_tree;
                const unsigned int              m_number_of_neighbors;
                std::vector<Element_index>      m_indices;
                Distance_adapter                m_distance_adapter;
            };

        } // namespace Region_growing_with_points 

    } // namespace Region_growing

} // namespace CGAL

#endif
