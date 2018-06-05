#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template <class Traits, class NormalMap>
            class Generalized_region_growing_points_conditions_2 {
            public:
                using Element       = Traits::Element;
                using Element_map   = Traits::Element_map;
                
                typedef NormalMap Normal_map;

                Generalized_region_growing_points_conditions_2(const Input_range& input_range) :
                    m_input_range(input_range)
                {}

                // Local condition
                bool is_in_same_region(const Element& elem1, const Element& elem2) const {
                    // ...
                    return true; // or false
                }

                // Global condition
                template <class Region>
                bool is_valid(const Region& region) const {
                    // return (region.size() > m_some_number);
                }
                
                // Update local conditions
                void update() {
                    // ...
                }
            private:
                const Input_range& m_input_range;
                const Element_map m_elem_map;
                const Normal_map m_normal_map;
            }
        }
    }
}
#endif
