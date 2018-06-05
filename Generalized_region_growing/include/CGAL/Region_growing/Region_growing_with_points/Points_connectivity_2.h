#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class Traits>
            class Generalized_region_growing_points_connectivity_2 {
            public:
                using Element     = Traits::Element;
                using Input_range = Traits::Input_range;
                using Kernel      = Traits::Kernel;
                using Element_map = Traits::Element_map;

                // ElementMap::key_type must be std::iterator_traits<Input_range::iterator>::value_type
                using Neighbors      = std::vector<Element_map::key_type>;
                using Neighbor_range = Iterator_range<Neighbors::const_iterator>;

                // Primary template
                template<class PointType, class ElementMapValueType>
                struct Search_structures {
                    using Search_base           = CGAL::Search_traits_2<Kernel>;
                    using Search_traits_adapter = CGAL::Search_traits_adapter<Element, Element_map, Search_base>;
                    using Fuzzy_circle          = CGAL::Fuzzy_sphere<Search_traits_adapter>;
                    using Tree                  = CGAL::Kd_tree<Search_traits_adapter>;
                }

                // Partial specialization when PointType and ElementMapValueType are the same
                template<class PointType>
                struct Search_structures<PointType, PointType> {
                    using Search_base   = CGAL::Search_traits_2<Kernel>;
                    using Fuzzy_circle  = CGAL::Fuzzy_sphere<Search_base>;
                    using Tree          = CGAL::Kd_tree<Search_base>;
                }

                using Point_2           = Kernel::Point_2;
                // Element_map::value_type must be Element
                using Search_structures = Search_structures<Point_2, Element_map::value_type>;
                using Fuzzy_circle      = Search_structures::Fuzzy_circle;
                using Tree              = Search_structures::Tree;

                Generalized_region_growing_points_connectivity_2(const Input_range& input_range) :
                    m_input_range(input_range),
                    m_tree(Tree(input_range.begin(), input_range.end()))
                { }

                void get_neighbors(const Element& center, Neighbor_range& output) const {
                    // The map must suit the input given
                    CGAL_precondition(std::is_same<
                                        Element_map::key_type, 
                                        std::iterator_traits<Input_range::iterator>::value_type
                                    >::value);
                    m_neighbors.clear();
                    Fuzzy_circle circle(center);
                    tree.search(std::back_inserter(m_neighbors), circle);
                    output = Neighbor_range(m_neighbors.begin(), m_neighbors.end());
                }
            private:
                Neighbors          m_neighbors;
                const Element_map  m_elem_map;
                const Input_range& m_input_range;
                Tree               m_tree;
            }
        }
    }
}

#endif
