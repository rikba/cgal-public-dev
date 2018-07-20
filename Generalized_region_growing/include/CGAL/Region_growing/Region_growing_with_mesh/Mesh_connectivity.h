#ifndef CGAL_REGION_GROWING_MESH_CONNECTIVITY_H
#define CGAL_REGION_GROWING_MESH_CONNECTIVITY_H

#include <Surface_mesh.h>
#include <Iterator_range.h>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            template <class Traits, class Mesh>
            class Mesh_connectivity {

            public:
                using Face = Mesh::Face_index;
                using Neighbor_range = CGAL::Iterator_range<CGAL::Face_around_face_iterator<Mesh> >;

                Meshes_connectivity(const Mesh& mesh) :
                m_mesh(mesh) { }

                Neighbors get_neighbors(const Face& face) {
                    return faces_around_face(m_mesh.halfedge(face), m_mesh);
                }

            private:
                const Mesh& m_mesh;
                Neighbor_range m_neighbor_range;
            };

        }
    }
}

#endif