#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H

template <class Traits>
class Generalized_region_growing_points_conditions_2 {
public:
    typedef typename Traits::Element Element;
    typedef typename Traits::Input_range::iterator Iterator;

    template <class ElementMap>
    bool is_same_region(Iterator iter1, Iterator iter2, const ElementMap& elem_map) const {
        Element elem1 = get(elem_map, iter1);
        Element elem2 = get(elem_map, iter2);
        // ...
        return true; // or false
    }
}

#endif