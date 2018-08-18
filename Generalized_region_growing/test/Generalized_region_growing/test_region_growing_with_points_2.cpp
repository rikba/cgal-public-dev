#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Timer.h>

std::ifstream in;

template < class Kernel >
bool test_region_growing_with_points_2(double radius = 2.9, size_t min_size = 5, double epsilon = 4.5, double normal_threshold = 0.7) {
    // default params for the data file inputbig_2.xyz

    // rewind to the beginning of the data file
    in.clear();
    in.seekg(0, std::ios::beg);

    using Point_2           = typename Kernel::Point_2;
    using Vector_2          = typename Kernel::Vector_2;
    using Point_with_normal = std::pair<Point_2, Vector_2>;
    using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
    using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
    using Input_range       = CGAL::Iterator_range<typename std::vector<Point_with_normal>::iterator>;
    using Traits            = CGAL::Region_growing::Region_growing_traits<Input_range, Point_map, Kernel>;
    using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_2<Traits, Normal_map>;
    using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_circular_query<Traits>;
    using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;
    using Region_range      = typename Region_growing::Region_range;
    using FT                = typename Kernel::FT;

    std::vector<Point_with_normal > pwn;
    double a, b, c, d, e, f;
    while (in >> a >> b >> c >> d >> e >> f) {
        pwn.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
    }

    Input_range input_range(pwn.begin(), pwn.end());
    Connectivity connectivity(input_range, radius);
    Conditions conditions(epsilon, normal_threshold, min_size);
    Region_growing region_growing(input_range, connectivity, conditions);
    region_growing.find_regions();
    Region_range regions = region_growing.regions();

    for (typename Region_range::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        if (!conditions.is_valid(*it)) return false;
    }
    return true;
}

int main(int argc, char *argv[]) {
    in = std::ifstream(argc > 1 ? argv[1] : "../data/inputbig_2.xyz");
    CGAL::set_ascii_mode(in);

    bool success = true;
    std::cout << "test_region_growing_with_points_2<CGAL::Exact_predicates_inexact_constructions_kernel>... ";
    if (!test_region_growing_with_points_2<CGAL::Exact_predicates_inexact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_2<CGAL::Exact_predicates_exact_constructions_kernel>... ";
    if (!test_region_growing_with_points_2<CGAL::Exact_predicates_exact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_2<CGAL::Simple_cartesian<float> >... ";
    if (!test_region_growing_with_points_2<CGAL::Simple_cartesian<float> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_2<CGAL::Simple_cartesian<double> >... ";
    if (!test_region_growing_with_points_2<CGAL::Simple_cartesian<double> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";
    
}
