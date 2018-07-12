#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/IO/write_ply_points.h>

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2           = Kernel::Point_2;
using Vector_2          = Kernel::Vector_2;
using Point_with_normal = std::pair<Point_2, Vector_2>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Input_range       = CGAL::Iterator_range<std::vector<Point_with_normal>::iterator>;

using Traits            = CGAL::Region_growing::Region_growing_with_points::Points_traits<Input_range, Point_map, Kernel>;
using Conditions        = CGAL::Region_growing::Region_growing_with_points::Points_conditions_2<Traits, Normal_map>;
using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity_circular_query<Traits>;
using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;

using Region_range      = Region_growing::Region_range;

using Point_3           = Kernel::Point_3;
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
            if (is_ascii(out))
                out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
            else
                out.write(reinterpret_cast<const char *>(&c), sizeof(c));
            return out;
        }
    };
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

    Connectivity connectivity(input_range, 2.9);

    Conditions conditions(input_range, 4.5, 0.7, 5);

    Region_growing region_growing(input_range, connectivity, conditions);

    clock_t start = clock();

    region_growing.find_regions();

    clock_t end = clock();

    Region_range regions = region_growing.regions();

    std::cerr << "Time elapsed: " << 1.0 * (end-start) / CLOCKS_PER_SEC << std::endl;

//    return 0;

    std::cout << std::setprecision(20);
    std::vector<Point_with_color> pwc;

    srand(time(NULL));
    for (Region_range::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        Color c = CGAL::make_array(static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256));
        const std::vector<int> &region = *it;
        for (int i : region) {
            Point_2 tmp = get(Point_map(), *(input_range.begin() + i));
            pwc.push_back(std::make_pair(Point_3(tmp.x(), tmp.y(), 0), c));
        }
    }
    CGAL::set_ascii_mode(std::cout);
    CGAL::write_ply_points_with_properties(std::cout, pwc,
                                           CGAL::make_ply_point_writer(PLY_Point_map()),
                                           std::make_tuple(PLY_Color_map(), CGAL::PLY_property<unsigned char>("red"),
                                                           CGAL::PLY_property<unsigned char>("green"),
                                                           CGAL::PLY_property<unsigned char>("blue")));

}
