#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class Traits, class Kernel>
            class Generalized_region_growing_points_connectivity_2 {
            public:
                typedef typename Traits::Element                         Element;
                typedef typename Traits::Input_range                     Input_range;
                typedef typename Traits::Element_map                     Element_map;

                typedef typename Input_range::const_iterator                Iterator;
                // TODO: Neighbors must be std::vector<Element> to adapt the kd tree's output
                typedef typename std::vector<Iterator>                      Neighbors;
                typedef typename Iterator_range<Neighbors::const_iterator>  Neighbor_range;

                // Primary template
                template<class PointType, class ElementMapValueType>
                struct Search_structures {
                    typedef CGAL::Search_traits_2<Kernel>                                  Search_base;
                    typedef CGAL::Search_traits_adapter<Element, Element_map, Search_base> Search_traits_adapter;
                    typedef CGAL::Fuzzy_sphere<Search_traits_adapter>                      Fuzzy_circle;
                    typedef CGAL::Kd_tree<Search_traits_adapter>                           Tree;
                }

                // Partial specialization when PointType and ElementMapValueType are the same
                template<class PointType>
                struct Search_structures<PointType, PointType> {
                    typedef CGAL::Search_traits_2<Kernel>   Search_base;
                    typedef CGAL::Fuzzy_sphere<Search_base> Fuzzy_circle;
                    typedef CGAL::Kd_tree<Search_base>      Tree;
                }

                typedef Kernel::Point_2                                     Point_2;
                typedef Search_structures<Point_2, Element_map::value_type> Search_structures;
                typedef Search_structures::Fuzzy_circle                     Fuzzy_circle;
                typedef Search_structures::Tree                             Tree;

                Generalized_region_growing_points_connectivity_2(const Input_range& input_range) :
                    m_input_range(input_range),
                    m_tree(Tree(input_range.begin(), input_range.end()))
                {}

                void get_neighbors(const Iterator& elem_iter, Neighbor_range& output) const {
                    m_neighbors.clear();
                    Element center = get(elem_map, *elem_iter);
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
