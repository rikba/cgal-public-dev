#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Region_growing.h>
#include <CGAL/property_map.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Surface_mesh.h>

using Kernel            = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3           = Kernel::Point_3;
using Vector_3          = Kernel::Vector_3;

using Mesh = CGAL::Surface_mesh<Point_3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;
using Faces = std::vector<Face_index>;
using Input_range = CGAL::Iterator_range<Faces::iterator>;
using Element_map = CGAL::Identity_property_map<Face_index>;

using Traits = CGAL::Region_growing::Region_growing_with_mesh::Mesh_traits<Input_range, Element_map, Kernel>;
using Conditions = CGAL::Region_growing::Region_growing_with_mesh::Mesh_conditions<Traits, Mesh>;
using Connectivity = CGAL::Region_growing::Region_growing_with_mesh::Mesh_connectivity<Traits, Mesh>;
using Region_growing = CGAL::Region_growing::Generalized_region_growing<Traits, Connectivity, Conditions>;


int main(int argc, char *argv[]) {
    Mesh m;
    Vertex_index v0 = m.add_vertex(Point_3(0, 0, 0));
    Vertex_index v1 = m.add_vertex(Point_3(0, 1, 0));
    Vertex_index v2 = m.add_vertex(Point_3(1, 0, 0));
    Vertex_index v3 = m.add_vertex(Point_3(1, 1, 0));
    Vertex_index v4 = m.add_vertex(Point_3(0, 0, 1));
    Vertex_index v5 = m.add_vertex(Point_3(0, 1, 1));
    Vertex_index v6 = m.add_vertex(Point_3(1, 0, 1));
    Vertex_index v7 = m.add_vertex(Point_3(1, 1, 1));
    m.add_face(v0, v3, v2);
    m.add_face(v0, v3, v1);
    m.add_face(v0, v6, v2);
    m.add_face(v0, v6, v4);
    m.add_face(v0, v5, v1);
    m.add_face(v0, v5, v4);
    m.add_face(v7, v2, v3);
    m.add_face(v7, v2, v6);
    m.add_face(v7, v1, v3);
    m.add_face(v7, v1, v5);
    m.add_face(v7, v4, v5);
    m.add_face(v7, v4, v6);

    Mesh::Face_range all_faces = m.faces();
    Faces faces;
    for (Mesh::Face_range::iterator it = all_faces.begin(); it != all_faces.end(); ++it)
        faces.push_back(*it);
    Input_range input_range(faces.begin(), faces.end());

    Connectivity connectivity(input_range, m);

    Conditions conditions(input_range, m, 0.1, 0.9);

    Region_growing rg(input_range, connectivity, conditions);

    rg.find_regions();

    std::cout << rg.number_of_regions() << std::endl;

}
