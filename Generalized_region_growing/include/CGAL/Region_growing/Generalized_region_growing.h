#ifndef CGAL_REGION_GROWING_GENERALIZED_REGION_GROWING_H
#define CGAL_REGION_GROWING_GENERALIZED_REGION_GROWING_H

#include <queue>
#include <CGAL/Iterator_range.h>

namespace CGAL {
    namespace Region_growing {

        template<class Traits, class Connectivity, class Conditions>
        class Generalized_region_growing {

        public:
            using Input_range             = typename Traits::Input_range;
            using Element                 = typename Traits::Element;
            using Element_map             = typename Traits::Element_map;
            using Element_with_properties = typename Element_map::key_type;

            using Neighbors               = std::list<Element_with_properties>;
            using Region                  = std::list<Element_with_properties>;
            using Regions                 = std::list<Region>;
            using Region_range            = CGAL::Iterator_range<typename Regions::const_iterator>;

            Generalized_region_growing(const Input_range &input_range, Connectivity &connectivity, Conditions &conditions) :
                m_input_range(input_range),
                m_connectivity(connectivity),
                m_conditions(conditions) { }

            void find_regions() {
                m_regions.clear();
                Region region;

                for (typename Input_range::const_iterator iter = m_input_range.begin(); iter != m_input_range.end(); ++iter) {

                    if (!m_visited[get(m_elem_map, *iter)]) { // Available element
                        region.clear();
                        // Grow a region from that element
                        grow_region(*iter, region);
                        // Check global condition
                        if (!m_conditions.is_valid(region)) {
                            revert(region);
                        } else {
                            m_regions.push_back(region);
                        }
                    }
                }

                m_output = Region_range(m_regions.begin(), m_regions.end());
            }

            const Region_range& regions() {
                return m_output;
            }

            const size_t number_of_regions() {
                return m_regions.size();
            }

        private:
            void grow_region(const Element_with_properties &seed, Region &region) {
                region.clear();
                // Use two queues, while running on this queue, push to the other queue;
                // When the queue is done, update the shape of the current region and swap to the other queue;
                // depth_index is the index of the queue we're using
                std::queue<Element_with_properties> ewp_queue[2];
                bool depth_index = 0;

                // Once an element is pushed to the queue, it is pushed to the region too.
                ewp_queue[depth_index].push(seed);
                m_visited[get(m_elem_map, seed)] = true;
                region.push_back(seed);

                m_conditions.update_shape(region);

                while (!ewp_queue[depth_index].empty() || !ewp_queue[!depth_index].empty()) {

                    // Call the next element of the queue and remove it from the queue
                    // but the element is not removed from the region.
                    Element_with_properties ewp = ewp_queue[depth_index].front();
                    ewp_queue[depth_index].pop();
                    // Get neighbors
                    Neighbors neighbors;
                    m_connectivity.get_neighbors(ewp, neighbors);

                    // Visit the neighbors;
                    for (typename Neighbors::const_iterator iter = neighbors.begin(); iter != neighbors.end(); iter++) {
                        Element_with_properties neighbor = *iter;

                        if (!m_visited[get(m_elem_map, neighbor)] && m_conditions.is_in_same_region(ewp, neighbor, region)) {

                            // Add to the other queue the neighbor which doesn't belong to any regions
                            // so that we can visit them later.
                            ewp_queue[!depth_index].push(neighbor);
                            region.push_back(neighbor);
                            m_visited[get(m_elem_map, neighbor)] = true;
                        }
                    }

                    if (ewp_queue[depth_index].empty()) {
                        m_conditions.update_shape(region);
                        depth_index = !depth_index;
                    }

                }
            }

            void revert(const Region& region) {
                for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it)
                    m_visited[get(m_elem_map, *it)] = false;
            }

        private:
            const Input_range &m_input_range;
            Element_map m_elem_map;
            Regions m_regions;
            Connectivity &m_connectivity;
            Conditions &m_conditions;
            std::map<Element, bool> m_visited;
            Region_range m_output = Region_range(m_regions.begin(), m_regions.end());
        };
    }
}

#endif
