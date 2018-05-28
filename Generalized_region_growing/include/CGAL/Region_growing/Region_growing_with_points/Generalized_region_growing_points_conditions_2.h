#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H

template <class Traits>
class Generalized_region_growing_points_conditions_2 {
public:
    typedef typename Traits::Element Element;
    typedef typename Traits::Element_map Element_map;
    typedef typename Traits::Normal_map Normal_map;
    typedef typename Traits::Input_range Input_range;
    typedef typename Input_range::iterator Iterator;

    Generalized_region_growing_points_conditions_2(Input_range input_range) :
        m_input_range(input_range)
    {}

    // Local condition
    bool is_in_same_region(Iterator elem_iter1, Iterator elem_iter2) const {
        Element elem1 = get(m_elem_map, elem_iter1);
        Element elem2 = get(m_elem_map, elem_iter2);
        // ...
        return true; // or false
    }

    // Global condition
    template <class Region>
    bool is_valid(const Region& region) const {
        // return (region.size() > m_some_number);
    }
private:
    Input_range m_input_range;
    Element_map m_elem_map;
    Normal_map m_normal_map;
}

#endif
