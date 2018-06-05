#ifndef CGAL_GENERALIZED_REGION_GROWING_H
#define CGAL_GENERALIZED_REGION_GROWING_H

namespace CGAL {
    namespace Region_growing {
        template <class Traits, class Connectivity, class Conditions>
        class Generalized_region_growing {
        public:
            using Element     = Traits::Element;
            using Input_range = Traits::Input_range;
            // Element_map::value_type is Element -- a geometry object
            using Element_map = Traits::Element_map;

            using Iterator      = InputRange::const_iterator;
            using Neighbors     = Connectivity::Neighbor_range;
            using Region        = std::vector<Element>;
            using Regions       = std::vector<Region>;
            using Region_range  = Iterator_range<Regions::const_iterator>;

            Generalized_region_growing(const Input_range& input_range, 
                                       const Connectivity& connectivity, 
                                       const Conditions& conditions) : 
                m_input_range(input_range),
                m_connectivity(connectivity),
                m_conditions(conditions)
            {}
            
            void find_regions() {
                for (Iterator iter = m_input_range.begin(); iter != m_input_range.end(); ++iter) {
                    if (!m_visited[get(m_elem_map, *iter)]) { // Available element
                        Region region;
                        // Grow a region from that element
                        grow_region(get(m_elem_map, *iter), region);
                        // Check global condition
                        if (m_conditions.is_valid(region))
                            m_regions.push_back(region); // Add the region grown
                        else {
                            // Revert the process
                            for (Region::const_iterator it = region.begin(); it != region.end(); ++it)
                                m_visited[*it] = false;
                        }
                    }
                }
                m_output = Region_range(m_regions.begin(), m_regions.end());
            }

            Region_range regions() {
                return m_output;
            }
        private:
            void grow_region(const Element& seed, Region& region) {
                region.clear();
                // Use a queue to keep track of the elements whose neighbors will be visited later.
                std::queue<Element> elem_queue;

                // Once an element is pushed to the queue, it is pushed to the region too.
                elem_queue.push(seed);
                m_visited[seed] = true;
                region.push(seed);

                while (!elem_queue.empty()) {

                    // Call the next element of the queue and remove it from the queue
                    // but the element is not removed from the region.
                    Element elem = elem_queue.front();
                    elem_queue.pop();

                    // Get neighbors
                    Neighbors neighbors;
                    m_connectivity.get_neighbors(elem, neighbors);

                    // Visit the neighbors;
                    for (Neighbors::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++) {
                        nb = get(m_elem_map, *iter); // neighbor
                        if (!m_visited[nb] && m_conditions.is_in_same_region(elem, nb)) {
                            // Add to the queue the neighbor which doesn't belong to any regions
                            // so that we can visit the neighbor's neighbors later.
                            elem_queue.push(nb);
                            region.push(nb);
                            m_visited[nb] = true;
                        }
                    }
                }
            }
        private:
            const Input_range& m_input_range;
            Regions m_regions;
            Region_range m_output;
            const Connectivity& m_connectivity;
            const Conditions& m_conditions;
            std::map<Element, bool> m_visited; // Keep track of available elements
        }
    }
}

#endif
