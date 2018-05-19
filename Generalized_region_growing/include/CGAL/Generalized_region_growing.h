#ifndef CGAL_GENERALIZED_REGION_GROWING_H
#define CGAL_GENERALIZED_REGION_GROWING_H

namespace CGAL {
    namespace Generalized_region_growing {
        template <class Traits, class Connectivity, class Conditions>
        class Generalized_region_growing {
        public:
            typedef typename Traits::Element Element;
            typedef typename Traits::InputRange Input_range;
            typedef typename InputRange::iterator Iterator;
            typedef typename Connectivity::OutputRange Neighbors;
            typedef typename std::vector<Iterator> Region;
            typedef typename std::vector<Region> Output;
            typedef typename std::map<Iterator, bool> State_map;

            template <class ElementMap>
            void find_regions(const Input_range& input_range, const ElementMap& elem_map) {
                for (Iterator iter = begin; iter != end; ++iter) {
                    if (!m_was_used[iter]) { // Available element
                        m_regions.push_back(Region()); // Add 1 region
                        Element elem = get(elem_map, iter);
                        // Grow a region from that element
                        grow_region(input_range, elem_map, m_state_map, iter, m_regions.back());
                    }
                }
            }
        private:
            template <class ElementMap>
            // Depth-first search from an element to form a region
            void grow_region(const Input_range& input_range, const ElementMap& elem_map, Iterator elem_iter, Region& region) {
                m_was_used[elem_iter] = true; // Mark the element as not available
                region.push_back(elem_iter); // Add the element to the region
                Neighbors neighbors;
                m_connectivity.get_neighbors(input_range, elem_map, elem_iter, neighbors); // Get neighbors
                for (Neighbors::iterator neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++) {
                    if (!m_was_used[iter] && m_conditions.is_same_region(elem_iter, neighbor_iter, elem_map)) {
                        // The neighbor is available and is in same region with the primary element
                        grow_region(input_range, elem_map, was_used, neighbor_iter, region);
                    }
                }
            }
        private:
            Output m_regions;
            Connectivity<Traits> m_connectivity;
            Conditions<Traits> m_conditions;
            State_map m_was_used; // Keep track of available elements
        }
    }
}

#endif