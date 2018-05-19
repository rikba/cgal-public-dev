#ifndef GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONNECTIVITY_2_H

#include <CGAL/Iterator_range.h>

template <class Traits>
class Generalized_region_growing_points_connectivity_2 {
public:
    typedef typename Traits::Element Element;
    typedef typename Traits::InputRange Input_range;
    typedef typename Input_range::iterator Iterator;
    typedef typename std::vector<Iterator> Output;
    typedef typename Iterator_range<Output::Iterator> OutputRange;

    template <class ElementMap>
    void get_neighbors(const Input_range& input, const ElementMap& elem_map, Iterator elem_iter, OutputRange& output) const {
        m_output.clear();
        Element elem1 = get(elem_map, elem_iter);
        /*
        for (Iterator iter = input.begin(); iter != input.end(); ++iter) {
            Element elem2 = get(elem_map, iter);
            ...
            m_output.push_back(iter);
        }
        */
       return OutputRange(m_output.begin(), m_output.end());
    }
private:
    Output m_output;
}

#endif
