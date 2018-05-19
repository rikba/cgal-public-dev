#ifndef GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_TRAITS_2_H

#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Iterator_range.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

class Generalized_region_growing_point_traits_2 {
public:
    typedef Kernel::Point_2 Element;
    typedef Iterator_range<std::vector<Element>::iterator> InputRange;
}

#endif