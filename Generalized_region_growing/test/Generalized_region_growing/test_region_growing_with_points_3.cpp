#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>

std::ifstream in;

template < class Kernel >
bool test_region_growing_with_points_3(size_t number_of_neighbors = 100, size_t min_size = 10, double epsilon = 0.5, double normal_threshold = 0.9) {
    // default params for the data file city_135.xyz

    // rewind to the beginning of the data file
    in.clear();
    in.seekg(0, std::ios::beg);

    using Point_3           = typename Kernel::Point_3;
    using Vector_3          = typename Kernel::Vector_3;
    using Point_with_normal = std::pair<Point_3, Vector_3>;
    using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
    using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
    using Input_range       = CGAL::Iterator_range<typename std::vector<Point_with_normal>::iterator>;
    using Traits            = CGAL::Region_growing::Region_growing_traits<Input_range, Point_map, Kernel>;
    using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_3<Traits, Normal_map>;
    using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_nearest_neighbors<Traits>;
    using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;
    using Region_range      = typename Region_growing::Region_range;

    std::vector<Point_with_normal> data;
    Point_3 point;
    Vector_3 normal;
    while (in >> point >> normal) {
        data.push_back(std::make_pair(point, normal));
    }

    Input_range input_range(data.begin(), data.end());
    Connectivity connectivity(input_range, number_of_neighbors);
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
    in = std::ifstream(argc > 1 ? argv[1] : "../data/city_135.xyz");
    CGAL::set_ascii_mode(in);

    bool success = true;
    std::cout << "test_region_growing_with_points_3<CGAL::Exact_predicates_inexact_constructions_kernel>... ";
    if (!test_region_growing_with_points_3<CGAL::Exact_predicates_inexact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_3<CGAL::Exact_predicates_exact_constructions_kernel>... ";
    if (!test_region_growing_with_points_3<CGAL::Exact_predicates_exact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_3<CGAL::Simple_cartesian<float> >... ";
    if (!test_region_growing_with_points_3<CGAL::Simple_cartesian<float> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_points_3<CGAL::Simple_cartesian<double> >... ";
    if (!test_region_growing_with_points_3<CGAL::Simple_cartesian<double> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";
    
}
