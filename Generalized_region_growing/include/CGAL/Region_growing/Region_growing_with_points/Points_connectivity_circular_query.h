#ifndef CGAL_GRG_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H
#define CGAL_GRG_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/property_map.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            template<class Traits>
            class Points_connectivity_circular_query {

            public:
                using Input_range             = typename Traits::Input_range;
                using Kernel                  = typename Traits::Kernel;
                using Element                 = typename Traits::Element; // Point_2 or Point_3
                using Element_map_original    = typename Traits::Element_map;
                
                struct Element_map {
                    using key_type = int;
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

                using Neighbors               = std::vector<int>;
                using Neighbor_range          = CGAL::Iterator_range<typename Neighbors::const_iterator>;

                using Search_base             = typename std::conditional<std::is_same<typename Kernel::Point_2, Element>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;
                using Search_traits_adapter   = CGAL::Search_traits_adapter<int, Element_map, Search_base>;
                using Fuzzy_circle            = CGAL::Fuzzy_sphere<Search_traits_adapter>;
                using Tree                    = CGAL::Kd_tree<Search_traits_adapter>;
                using Splitter                = CGAL::Sliding_midpoint<Search_traits_adapter>;

                using FT                      = typename Kernel::FT;

                Points_connectivity_circular_query(const Input_range &input_range, const FT &radius) :
                m_input_range(input_range),
                m_radius(radius),
                m_elem_map(Element_map(input_range)),
                m_sta(Search_traits_adapter(Element_map(input_range), Search_base())),
                m_tree(Splitter(), Search_traits_adapter(Element_map(input_range), Search_base())) {

                    for(int i = 0; i < m_input_range.end() - m_input_range.begin(); ++i)
                        indices.push_back(i);

                    m_tree.insert(indices.begin(), indices.end());

                }

                void get_neighbors(const int &center, Neighbors& neighbors) {

                    neighbors.clear();
                    Fuzzy_circle circle(center, m_radius, 0, m_sta);
                    m_tree.search(std::back_inserter(neighbors), circle);

                }

            private:
                const Element_map               m_elem_map;
                const Input_range &             m_input_range;
                Tree                            m_tree;
                const FT                        m_radius;
                const Search_traits_adapter     m_sta;
                Neighbors                       indices;
            };

        } // namespace Region_growing_with_points 

    } // namespace Region_growing

} // namespace CGAL

#endif
