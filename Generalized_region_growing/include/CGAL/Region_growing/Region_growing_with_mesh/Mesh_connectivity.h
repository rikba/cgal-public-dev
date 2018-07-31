#ifndef CGAL_REGION_GROWING_MESH_CONNECTIVITY_H
#define CGAL_REGION_GROWING_MESH_CONNECTIVITY_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Iterator_range.h>
#include <sstream>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            template <class Traits, class Mesh>
            class Mesh_connectivity {

            public:
                using Face                    = typename Traits::Element; // Mesh::Face_index
                using Input_range             = typename Traits::Input_range;
                using Face_range              = CGAL::Iterator_range<CGAL::Face_around_face_iterator<Mesh> >;
                using Element_map             = typename Traits::Element_map;
                using Element_with_properties = typename Element_map::key_type;

                Mesh_connectivity(const Mesh& mesh) :
                m_mesh(mesh) { }

                template < class Neighbors_ >
                void get_neighbors(const Element_with_properties& ewp, Neighbors_& neighbors) {
                    Face face = get(m_elem_map, ewp);
                    neighbors.clear();
                    const Face_range& tmp = faces_around_face(m_mesh.halfedge(face), m_mesh);
                    for (typename Face_range::iterator it = tmp.begin(); it != tmp.end(); ++it) {
                        neighbors.push_back(*it);
                    }
                }

            private:
                const Mesh& m_mesh;
                const Element_map m_elem_map = Element_map();
            };

        }
    }
}

#endif