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
                using Face                  = typename Mesh::Face_index;
                using Input_range           = typename Traits::Input_range;
                using Face_range            = CGAL::Iterator_range<CGAL::Face_around_face_iterator<Mesh> >;
                using Element_index         = size_t;

                Mesh_connectivity(const Input_range& input_range, const Mesh& mesh) :
                m_mesh(mesh),
                m_input_range(input_range) { }

                template < class Neighbors_ >
                void get_neighbors(Element_index element_index, Neighbors_& neighbors) {
                    Face face = *(m_input_range.begin() + element_index);
                    neighbors.clear();
                    const Face_range& tmp = faces_around_face(m_mesh.halfedge(face), m_mesh);
                    for (typename Face_range::iterator it = tmp.begin(); it != tmp.end(); ++it) {
                        std::stringstream stream;
                        stream.clear();
                        stream.str("");
                        stream << (*it);
                        char ch;
                        size_t nb;
                        stream >> ch >> nb;
                        if (nb != 4294967295) neighbors.push_back(nb);
                    }
                }

            private:
                const Mesh& m_mesh;
                const Input_range& m_input_range;
            };

        }
    }
}

#endif