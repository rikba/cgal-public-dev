#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Timer.h>

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

using Point_3           = Kernel::Point_3;
using Color             = CGAL::cpp11::array<unsigned char, 3>;
using Point_with_color  = std::pair<Point_3, Color>;
using PLY_Point_map     = CGAL::First_of_pair_property_map<Point_with_color>;
using PLY_Color_map     = CGAL::Second_of_pair_property_map<Point_with_color>;

using FT                = Kernel::FT;

void benchmark_region_growing(const size_t test_count, const Input_range& input_range, const double radius, const size_t min_size, const double epsilon, const double normal_threshold) {

    Connectivity connectivity(input_range, radius);
    Conditions conditions(epsilon, normal_threshold, min_size);
    Region_growing region_growing(input_range, connectivity, conditions);

    CGAL::Timer t;
    t.start();
    region_growing.find_regions();
    t.stop();

    size_t number_of_points_assigned = 0;
    Region_growing::Region_range regions = region_growing.regions();

    for (Region_range::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        number_of_points_assigned += (*it).size();
    }

    std::cout << "Test #" << test_count << ":\n";
    std::cout << "  radius = " << radius << ";\n";
    std::cout << "  min_size = " << min_size << ";\n";
    std::cout << "  epsilon = " << epsilon << ";\n";
    std::cout << "  normal_threshold = " << normal_threshold << ";\n";
    std::cout << "  -----\n";
    std::cout << "  Time elapsed: " << t.time() << '\n';
    std::cout << "  Number of regions detected: " << region_growing.number_of_regions() << '\n';
    std::cout << "  Number of points assigned: " << number_of_points_assigned << '\n';
    std::cout << '\n';
}

int main(int argc, char *argv[]) {
    std::ifstream in(argc > 1 ? argv[1] : "../data/inputbig_2.xyz");
    CGAL::set_ascii_mode(in);

    std::vector<Point_with_normal > pwn;
    double a,b,c,d,e,f;
    int i = 0;
    while (in >> a >> b >> c >> d >> e >> f) {
        pwn.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
        ++i;
    }

    Input_range input_range(pwn.begin(), pwn.end());

    benchmark_region_growing(1, input_range, 1, 5, 4.5, 0.7);
    benchmark_region_growing(2, input_range, 3, 5, 4.5, 0.7);
    benchmark_region_growing(3, input_range, 6, 5, 4.5, 0.7);
    benchmark_region_growing(4, input_range, 9, 5, 4.5, 0.7);
}
