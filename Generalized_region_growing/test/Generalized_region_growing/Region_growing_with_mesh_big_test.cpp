#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3           = Kernel::Point_3;
using Vector_3          = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;
using Faces = std::vector<Face_index>;
using Input_range = CGAL::Iterator_range<Faces::iterator>;
using Element_map = CGAL::Identity_property_map<Face_index>;

using Traits = CGAL::Region_growing::Region_growing_traits<Input_range, Element_map, Kernel>;
using Conditions = CGAL::Region_growing::Region_growing_with_mesh::Mesh_conditions<Traits, Mesh>;
using Connectivity = CGAL::Region_growing::Region_growing_with_mesh::Mesh_connectivity<Traits, Mesh>;
using Region_growing = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;


int main(int argc, char *argv[]) {
    Mesh m;
    std::ifstream in(argc > 1 ? argv[1] : "../data/car_mesh.inp");

    int vertex_count, face_count;
    in >> vertex_count;
    std::vector<Mesh::Vertex_index> vv;
    for (int i = 0; i < vertex_count; ++i) {
        double x, y, z;
        in >> x >> y >> z;
        vv.push_back(m.add_vertex(Point_3(x, y, z)));
    }
    in >> face_count;
    Faces faces;
    std::vector<boost::tuple<int, int, int> > ply_faces;
    for (int i = 0; i < face_count; ++i) {
        int v1, v2, v3;
        in >> v1 >> v2 >> v3;
        if (CGAL::collinear(m.point(vv[v1-1]), m.point(vv[v2-1]), m.point(vv[v3-1]))) continue;
        faces.push_back(m.add_face(vv[v1-1], vv[v2-1], vv[v3-1]));
        ply_faces.push_back(boost::make_tuple(v1-1, v2-1, v3-1));
    }

    Input_range input_range(faces.begin(), faces.end());

    Connectivity connectivity(m);

    double epsilon, normal_threshold;
    std::cerr << "Input epsilon: ";
    std::cin >> epsilon;
    std::cerr << "Input normal_threshold: ";
    std::cin >> normal_threshold;

    // epsilon = 5, normal_threshold = 0.7

    Conditions conditions(m, epsilon, normal_threshold);

    Region_growing rg(input_range, connectivity, conditions);

    std::vector<size_t> v;

    CGAL::Timer t;
    t.start();

    rg.find_regions();

    t.stop();

    std::cerr << "comment Time elapsed: " << 1.0 * clock() / CLOCKS_PER_SEC << " s.\n";
    std::cerr << "comment " << rg.number_of_regions() << " regions found." << '\n';

    std::cout << "ply\nformat ascii 1.0\n";

    std::vector<boost::tuple<Point_3, int, int, int> > ply_vertices;

    Region_growing::Region_range all_regions = rg.regions();

    srand(time(NULL));

    for (Region_growing::Region_range::iterator it = all_regions.begin(); it != all_regions.end(); ++it) {
        int r = rand()%256, g = rand()%256, b = rand()%256;
        if ((*it).size() < 5) r = g = b = 255;
        for (size_t face_id : *it) {
            auto vertices = m.vertices_around_face(m.halfedge(faces[face_id]));
            for (auto itv = vertices.begin(); itv != vertices.end(); ++itv) {
                ply_vertices.push_back(boost::make_tuple(m.point(*itv), r, g, b));
            }
            int sz = (int)ply_vertices.size();
            ply_faces[face_id] = boost::make_tuple(sz-1, sz-2, sz-3);
        }
    }

    std::cout << "element vertex " << ply_vertices.size() << "\n" << "property float x\nproperty float y\nproperty float z\n";
    std::cout << "property uchar red\nproperty uchar green\nproperty uchar blue\n";

    std::cout << "element face " << ply_faces.size() << "\nproperty list uchar int vertex_index\nend_header\n";
    for (int i = 0; i < ply_vertices.size(); ++i)
        std::cout << boost::get<0>(ply_vertices[i]) << ' ' << boost::get<1>(ply_vertices[i]) << ' ' << boost::get<2>(ply_vertices[i]) << ' ' << boost::get<3>(ply_vertices[i]) << '\n';
//    return 0;
    for (int i = 0; i < ply_faces.size(); ++i)
        std::cout << 3 << ' ' << boost::get<0>(ply_faces[i]) << ' ' << boost::get<1>(ply_faces[i]) << ' ' << boost::get<2>(ply_faces[i]) << '\n';

}
