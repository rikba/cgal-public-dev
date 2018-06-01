#ifndef CGAL_GENERALIZED_REGION_GROWING_H
#define CGAL_GENERALIZED_REGION_GROWING_H

namespace CGAL {
    namespace Region_growing {
        template <class Traits, class Connectivity, class Conditions>
        class Generalized_region_growing {
        public:
            typedef typename Traits::Element Element;
            typedef typename Traits::InputRange Input_range;
            typedef typename Traits::Element_map Element_map;

            typedef typename InputRange::const_iterator Iterator;
            typedef typename Connectivity::Neighbor_range Neighbors;
            typedef typename std::vector<Iterator> Region;
            typedef typename std::vector<Region> Regions;
            typedef typename Iterator_range<Regions::const_iterator> Region_range;

            Generalized_region_growing(const Input_range& input_range, 
                                        const Connectivity& connectivity, 
                                        const Conditions& conditions) : 
                m_input_range(input_range),
                m_connectivity(connectivity),
                m_conditions(conditions)
            {}
            
            void find_regions() {
                for (Iterator iter = m_input_range.begin(); iter != m_input_range.end(); ++iter) {
                    if (!m_visited[iter]) { // Available element
                        Region region;
                        Element elem = get(m_elem_map, *iter);
                        // Grow a region from that element
                        grow_region(iter, region);
                        // Check global condition
                        if (m_conditions.is_valid(region))
                            m_regions.push_back(region); // Add the region grown
                        else {
                            // Revert the process
                            // NOTE: Region::const_iterator is a double iterator
                            for (Region::const_iterator it = region.begin(); it != region.end(); ++it)
                                m_visited[*(*it)] = false;
                        }
                    }
                }
                m_output = Region_range(m_regions.begin(), m_regions.end());
            }

            Region_range regions() {
                return m_output;
            }
        private:
            void grow_region(const Iterator& seed_iter, Region& region) {
                region.clear();
                // Use a queue to keep track of the elements whose neighbors will be visited later.
                std::queue<Iterator> elem_queue;

                // Once an element is pushed to the queue, it is pushed to the region too.
                elem_queue.push(seed_iter);
                m_visited[elem_iter] = true;
                region.push(elem_iter);

                while (!elem_queue.empty()) {

                    // Call the next element of the queue and remove it from the queue
                    // but the element is not removed from the region.
                    Iterator elem_iter = elem_queue.front();
                    elem_queue.pop();

                    // Get neighbors
                    Neighbors neighbors;
                    m_connectivity.get_neighbors(elem_iter, neighbors);

                    // Visit the neighbors;
                    // NOTE: Neighbors::const_iterator is a double pointer
                    for (Neighbors::const_iterator nb_iter = neighbors.begin(); nb_iter != neighbors.end(); nb_iter++) {
                        if (!m_visited[*nb_iter] && m_conditions.is_in_same_region(elem_iter, *nb_iter, elem_map)) {
                            // Add to the queue the neighbor which doesn't belong to any regions
                            // so that we can visit the neighbor's neighbors later.
                            elem_queue.push(*nb_iter);
                            region.push(*nb_iter);
                            m_visited[*nb_iter] = true;
                        }
                    }
                }
            }
        private:
            const Input_range& m_input_range;
            const Element_map m_elem_map;
            Regions m_regions;
            Region_range m_output;
            const Connectivity& m_connectivity;
            const Conditions& m_conditions;
            std::map<Iterator, bool> m_visited; // Keep track of available elements
        }
    }
}

#endif
