#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>

// Type declaration

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2           = Kernel::Point_2;
using Vector_2          = Kernel::Vector_2;
using Point_with_normal = std::pair<Point_2, Vector_2>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Input_range       = CGAL::Iterator_range<std::vector<Point_with_normal>::iterator>;

using Traits            = CGAL::Region_growing::Region_growing_traits<Input_range, Point_map, Kernel>;
using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_2<Traits, Normal_map>;
using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_circular_query<Traits>;
using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;

using Region_range      = Region_growing::Region_range;

int main(int argc, char *argv[]) {
    // Prepare the input
    std::ifstream in(argc > 1 ? argv[1] : "../data/inputbig_2.xyz");
    CGAL::set_ascii_mode(in);
    std::vector<Point_with_normal > pwn;
    double a,b,c,d,e,f;
    while (in >> a >> b >> c >> d >> e >> f) {
        pwn.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
    }
    Input_range input_range(pwn.begin(), pwn.end());

    // Create instances of the class Connectivity and Conditions
    Connectivity connectivity(input_range, 5);
    Conditions conditions(4.5, 0.7, 5);

    // Create an instance of the region growing class
    Region_growing region_growing(input_range, connectivity, conditions);

    // Run the algorithm
    region_growing.find_regions();

    // Print the number of regions found
    std::cerr << region_growing.number_of_regions() << " regions found." << '\n';

}
