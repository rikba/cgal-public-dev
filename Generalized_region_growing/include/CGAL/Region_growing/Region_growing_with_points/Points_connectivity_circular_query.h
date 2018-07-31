#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_CIRCULAR_QUERY_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            template<class Traits>
            class Points_connectivity_circular_query {

            public:
                using Input_range             = typename Traits::Input_range;
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;
                using Element                 = typename Traits::Element;
                using Element_with_properties = typename Element_map::key_type;

                using Search_base             = typename std::conditional<std::is_same<typename Kernel::Point_2, Element>::value, CGAL::Search_traits_2<Kernel>, CGAL::Search_traits_3<Kernel> >::type;
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

                // Element == Point_3 or Point_2
                using Point             = typename Traits::Element;
                using FT                = typename Kernel::FT;
                using Kd_tree_config    = Search_structures<Point, Element_with_properties>;
                using Fuzzy_circle      = typename Kd_tree_config::Fuzzy_circle;
                using Tree              = typename Kd_tree_config::Tree;

                Points_connectivity_circular_query(const Input_range &input_range, const FT radius) :
                m_input_range(input_range),
                m_radius(radius) {
                    m_tree.insert(m_input_range.begin(), m_input_range.end());
                }

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