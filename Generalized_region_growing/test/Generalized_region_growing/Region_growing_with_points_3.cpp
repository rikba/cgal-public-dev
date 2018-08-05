#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/random_simplify_point_set.h>

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3           = Kernel::Point_3;
using Vector_3          = Kernel::Vector_3;
using Point_with_normal = std::pair<Point_3, Vector_3>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Input_range       = CGAL::Iterator_range<std::vector<Point_with_normal>::iterator>;

using Traits            = CGAL::Region_growing::Region_growing_traits<Input_range, Point_map, Kernel>;
using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_3<Traits, Normal_map>;
using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_nearest_neighbors<Traits>;
using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;

using Region_range      = Region_growing::Region_range;

using Color             = CGAL::cpp11::array<unsigned char, 3>;
using Point_with_color  = std::pair<Point_3, Color>;
using PLY_Point_map     = CGAL::First_of_pair_property_map<Point_with_color>;
using PLY_Color_map     = CGAL::Second_of_pair_property_map<Point_with_color>;

// Define how a color should be stored
namespace CGAL {
    template<class F>
    struct Output_rep<::Color, F> {
        const ::Color &c;
        static const bool is_specialized = true;

        Output_rep(const ::Color &c) : c(c) {}

        std::ostream &operator()(std::ostream &out) const {
            if (CGAL::is_ascii(out))
                out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
            else
                out.write(reinterpret_cast<const char *>(&c), sizeof(c));
            return out;
        }
    };
}

int main(int argc, char *argv[]) {
    std::ifstream in(argc > 1 ? argv[1] : "../data/city_135.xyz");
    CGAL::set_ascii_mode(in);

    std::vector<Point_with_normal> data;
    Point_3 point;
    Vector_3 normal;
    while (in >> point >> normal) {
        data.push_back(std::make_pair(point, normal));
    }

    Input_range input_range(data.begin(), data.end());

    int number_of_neighbors, min_size;
    double epsilon, normal_threshold;

    std::cerr << "Input number_of_neighbors: "; std::cin >> number_of_neighbors;
    std::cerr << "Input min_size: "; std::cin >> min_size;
    std::cerr << "Input epsilon: "; std::cin >> epsilon;
    std::cerr << "Input normal_threshold: "; std::cin >> normal_threshold;

    // number_of_neighbors = 100;
    // min_size = 1;
    // epsilon = 0.5;
    // normal_threshold = 0.9;

    Connectivity connectivity(input_range, number_of_neighbors);

    Conditions conditions(epsilon, normal_threshold, min_size);

    Region_growing region_growing(input_range, connectivity, conditions);

    region_growing.find_regions();

    Region_range regions = region_growing.regions();

    std::cerr << "comment Time elapsed: " << 1.0 * clock() / CLOCKS_PER_SEC << " s.\n";
    std::cerr << "comment " << region_growing.number_of_regions() << " regions found." << '\n';

    std::vector<Point_with_color> pwc;

    srand(time(NULL));
    for (Region_range::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        Color c = CGAL::make_array(static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256));
        const Region_growing::Region &region = *it;
        for (Region_growing::Element_with_properties ewp : region) {
            Point_3 tmp = get(Point_map(), ewp);
            pwc.push_back(std::make_pair(tmp, c));
        }
    }

    std::cerr << "comment " << pwc.size() << " points assigned." << '\n';

    CGAL::set_ascii_mode(std::cout);
    CGAL::write_ply_points_with_properties(std::cout, pwc,
                                           CGAL::make_ply_point_writer(PLY_Point_map()),
                                           std::make_tuple(PLY_Color_map(), CGAL::PLY_property<unsigned char>("red"),
                                                           CGAL::PLY_property<unsigned char>("green"),
                                                           CGAL::PLY_property<unsigned char>("blue")));

}
