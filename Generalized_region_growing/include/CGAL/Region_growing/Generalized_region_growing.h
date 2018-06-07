#include <queue>
#include <CGAL/Iterator_range.h>

#ifndef CGAL_GENERALIZED_REGION_GROWING_H
#define CGAL_GENERALIZED_REGION_GROWING_H

namespace CGAL {
    namespace Region_growing {
        template<class Traits, class Connectivity, class Conditions>
        class Generalized_region_growing {
        public:
            using Input_range             = typename Traits::Input_range;
            using Element                 = typename Traits::Element;
            using Element_map             = typename Traits::Element_map;
            using Element_with_properties = typename Element_map::key_type;

            using Iterator      = typename Input_range::const_iterator;
            using Neighbors     = typename Connectivity::Neighbor_range;
            using Region        = std::vector<Element>;
            using Regions       = std::vector<Region>;
            using Region_range  = CGAL::Iterator_range<typename Regions::const_iterator>;

            Generalized_region_growing(const Input_range &input_range,
                                       Connectivity &connectivity,
                                       Conditions &conditions) :
                    m_input_range(input_range),
                    m_connectivity(connectivity),
                    m_conditions(conditions) {}

            void find_regions() {
                for (Iterator iter = m_input_range.begin(); iter != m_input_range.end(); ++iter) {
                    Element_with_properties ewp = *iter;
                    Element elem = m_elem_map[ewp];
                    if (!m_visited[elem]) { // Available element
                        Region region;
                        // Grow a region from that element
                        grow_region(ewp, region);
                        // Check global condition
                        if (m_conditions.is_valid(region))
                            m_regions.push_back(region); // Add the region grown
                        else {
                            // Revert the process
                            for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it)
                                m_visited[*it] = false;
                        }
                    }
                }
            }

            Region_range regions() {
                return Region_range(m_regions.begin(), m_regions.end());
            }

            void print() {
                for (Region r : m_regions) {
                    for (Element e : r) {
                        std::cout << e << " ";
                    }
                    std::cout << '\n';
                }
            }

        private:
            void grow_region(const Element_with_properties &seed_with_properties, Region &region) {
                region.clear();
                // Use a queue to keep track of the elements whose neighbors will be visited later.
                // The queue keeps properties too, but region and state map don't use it.
                std::queue<Element_with_properties> ewp_queue;
                Element seed = get(m_elem_map, seed_with_properties);

                // Once an element is pushed to the queue, it is pushed to the region too.
                ewp_queue.push(seed_with_properties);
                m_visited[seed] = true;
                region.push_back(seed);

                while (!ewp_queue.empty()) {

                    // Call the next element of the queue and remove it from the queue
                    // but the element is not removed from the region.
                    Element_with_properties ewp = ewp_queue.front();
                    ewp_queue.pop();

                    // Get neighbors
                    Neighbors neighbors = m_connectivity.get_neighbors(ewp);

                    // Visit the neighbors;
                    for (typename Neighbors::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++) {
                        Element_with_properties nwp = *iter; // neighbor with properties
                        Element neighbor = m_elem_map[nwp]; // neighbor

                        if (!m_visited[neighbor] && m_conditions.is_in_same_region(ewp, nwp, region)) {
                            // Add to the queue the neighbor which doesn't belong to any regions
                            // so that we can visit the neighbor's neighbors later.
                            ewp_queue.push(nwp);
                            region.push_back(neighbor);
                            m_visited[neighbor] = true;
                        }
                    }
                }
            }

        private:
            const Input_range &m_input_range;
            Element_map m_elem_map;
            Regions m_regions;
            Connectivity &m_connectivity;
            Conditions &m_conditions;
            std::map<Element, bool> m_visited; // Keep track of available elements
        };
    }
}

#endif
