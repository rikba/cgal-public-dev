#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
using Connectivity      = CGAL::Region_growing::Region_growing_with_points::Points_connectivity<Traits>;
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
    std::ifstream in(argc > 1 ? argv[1] : "../data/quadrilateral_and_circle.xy");
    CGAL::set_ascii_mode(in);

    Point_2 point;
    Vector_2 normal;

    std::vector<Point_with_normal> pwn;
    while (in >> point >> normal) {
        pwn.push_back(std::make_pair(point, normal));
    }

    Input_range input_range(pwn.begin(), pwn.end());

    Connectivity connectivity(input_range, 0.5);

    Conditions conditions(0.2, 0.3, 0);

    Region_growing region_growing(input_range, connectivity, conditions);

    region_growing.find_regions();

    Region_range regions = region_growing.regions();

    std::vector<Point_with_color> pwc;

    srand(time(NULL));
    for (Region_range::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        Color c = CGAL::make_array(static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256),
                                   static_cast<unsigned char>(rand() % 256));
        const std::vector<Point_with_normal> &region = *it;
        for (Point_with_normal p : region) {
            pwc.push_back(std::make_pair(Point_3(p.first.x(), p.first.y(), 0), c));
        }
    }
    CGAL::set_ascii_mode(std::cout);
    CGAL::write_ply_points_with_properties(std::cout, pwc,
                                           CGAL::make_ply_point_writer(PLY_Point_map()),
                                           std::make_tuple(PLY_Color_map(), CGAL::PLY_property<unsigned char>("red"),
                                                           CGAL::PLY_property<unsigned char>("green"),
                                                           CGAL::PLY_property<unsigned char>("blue")));

}
