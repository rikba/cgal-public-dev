#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

std::ifstream in;

template < class Kernel >
bool test_region_growing_with_mesh(double epsilon = 5.0, double normal_threshold = 0.7) { 
    // default param value for car_mesh.inp

    // rewind to the beginning of the data file
    in.clear();
    in.seekg(0, std::ios::beg);

    using Point_3           = typename Kernel::Point_3;
    using Vector_3          = typename Kernel::Vector_3;
    using Mesh              = CGAL::Surface_mesh<Point_3>;
    using Vertex_index      = typename Mesh::Vertex_index;
    using Face_index        = typename Mesh::Face_index;
    using Faces             = std::vector<Face_index>;
    using Input_range       = CGAL::Iterator_range<typename Faces::iterator>;
    using Element_map       = CGAL::Identity_property_map<Face_index>;
    using Traits            = CGAL::Region_growing::Region_growing_traits<Input_range, Element_map, Kernel>;
    using Conditions        = CGAL::Region_growing::Region_growing_with_mesh::Mesh_conditions<Traits, Mesh>;
    using Connectivity      = CGAL::Region_growing::Region_growing_with_mesh::Mesh_connectivity<Traits, Mesh>;
    using Region_growing    = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;

    Mesh m;

    int vertex_count, face_count;
    in >> vertex_count;
    std::vector<typename Mesh::Vertex_index> vv;
    for (int i = 0; i < vertex_count; ++i) {
        double x, y, z;
        in >> x >> y >> z;
        vv.push_back(m.add_vertex(Point_3(x, y, z)));
    }
    in >> face_count;
    Faces faces;
    for (int i = 0; i < face_count; ++i) {
        int v1, v2, v3;
        in >> v1 >> v2 >> v3;
        if (CGAL::collinear(m.point(vv[v1-1]), m.point(vv[v2-1]), m.point(vv[v3-1]))) continue;
        faces.push_back(m.add_face(vv[v1-1], vv[v2-1], vv[v3-1]));
    }

    Input_range input_range(faces.begin(), faces.end());
    Connectivity connectivity(m);
    Conditions conditions(m, epsilon, normal_threshold);

    Region_growing rg(input_range, connectivity, conditions);
    rg.find_regions();
    typename Region_growing::Region_range all_regions = rg.regions();

    for (typename Region_growing::Region_range::iterator it = all_regions.begin(); it != all_regions.end(); ++it) {
        if (!conditions.is_valid(*it)) return false;
    }
    return true;
}


int main(int argc, char *argv[]) {
    in = std::ifstream(argc > 1 ? argv[1] : "../data/car_mesh.inp");
    bool success = true;
    std::cout << "test_region_growing_with_mesh<CGAL::Exact_predicates_inexact_constructions_kernel>... ";
    if (!test_region_growing_with_mesh<CGAL::Exact_predicates_inexact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_mesh<CGAL::Exact_predicates_exact_constructions_kernel>... ";
    if (!test_region_growing_with_mesh<CGAL::Exact_predicates_exact_constructions_kernel>()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_mesh<CGAL::Simple_cartesian<float> >... ";
    if (!test_region_growing_with_mesh<CGAL::Simple_cartesian<float> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";

    std::cout << "test_region_growing_with_mesh<CGAL::Simple_cartesian<double> >... ";
    if (!test_region_growing_with_mesh<CGAL::Simple_cartesian<double> >()) {
        success = false;
        std::cout << "failed\n";
    } else std::cout << "succeeded\n";
}
