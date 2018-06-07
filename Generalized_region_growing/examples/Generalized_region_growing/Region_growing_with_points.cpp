#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2           = Kernel::Point_2;
using Vector_2          = Kernel::Vector_2;
using Point_with_normal = std::pair<Point_2, Vector_2>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Input_range       = CGAL::Iterator_range<std::vector<Point_with_normal>::iterator>;

using Traits            = CGAL::Region_growing::Region_growing_with_points::Points_traits_2<Input_range, Point_map, Kernel>;
using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_2<Traits, Normal_map>;
using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_2<Traits>;
using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;

int main() {
    std::ifstream stream("../data/points.xy");
    CGAL::set_ascii_mode(stream);

    Point_2 point;
    Vector_2 normal;

    std::vector<Point_with_normal> pwn;
    while (stream >> point >> normal) {
        pwn.push_back(std::make_pair(point, normal));
    }

    Input_range input_range(pwn.begin(), pwn.end());

    Connectivity connectivity(input_range, 2.0);

    Conditions conditions(1.0, 0.5, 0);

    Region_growing region_growing(input_range, connectivity, conditions);

    region_growing.find_regions();

    region_growing.print();

}
