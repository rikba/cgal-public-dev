#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template<class Traits>
            class Points_connectivity_2 {
            public:
                using Element                 = typename Traits::Element;
                using Input_range             = typename Traits::Input_range;
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;

                using Element_with_properties = typename Element_map::key_type;

                using Neighbors      = std::vector<Element_with_properties>;
                using Neighbor_range = CGAL::Iterator_range<typename Neighbors::const_iterator>;

                // Primary template
                template<class PointType, class ElementWithProperties>
                struct Search_structures {
                    using Search_base           = CGAL::Search_traits_2<Kernel>;
                    using Search_traits_adapter = CGAL::Search_traits_adapter<Element_with_properties, Element_map, Search_base>;
                    using Fuzzy_circle          = CGAL::Fuzzy_sphere<Search_traits_adapter>;
                    using Tree                  = CGAL::Kd_tree<Search_traits_adapter>;
                };

                // Partial specialization when PointType and ElementWithProperties are the same
                template<class PointType>
                struct Search_structures<PointType, PointType> {
                    using Search_base   = CGAL::Search_traits_2<Kernel>;
                    using Fuzzy_circle  = CGAL::Fuzzy_sphere<Search_base>;
                    using Tree          = CGAL::Kd_tree<Search_base>;
                };

                // Element == Point_2
                using Point_2           = typename Kernel::Point_2;
                using FT                = typename Kernel::FT;
                using Kd_tree_config    = Search_structures<Point_2, Element_with_properties>;
                using Fuzzy_circle      = typename Kd_tree_config::Fuzzy_circle;
                using Tree              = typename Kd_tree_config::Tree;

                Points_connectivity_2(const Input_range &input_range, const FT &radius) :
                        m_input_range(input_range),
                        m_radius(radius) {
                    m_tree.insert(m_input_range.begin(), m_input_range.end());
                }

                Neighbor_range get_neighbors(const Element_with_properties &center) {
//                    Tree tree(m_input_range.begin(), m_input_range.end());
                    m_neighbors.clear();
                    Fuzzy_circle circle(center, m_radius);
                    m_tree.search(std::back_inserter(m_neighbors), circle);
                    return Neighbor_range(m_neighbors.begin(), m_neighbors.end());
                }

            private:
                Neighbors m_neighbors;
                const Element_map m_elem_map = Element_map();
                const Input_range &m_input_range;
                Tree m_tree;
                const FT m_radius;
            };
        }
    }
}

#endif
