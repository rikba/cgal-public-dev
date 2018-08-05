#ifndef CGAL_REGION_GROWING_MESH_CONNECTIVITY_H
#define CGAL_REGION_GROWING_MESH_CONNECTIVITY_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Iterator_range.h>
#include <sstream>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            /*!
                \ingroup PkgGeneralizedRegionGrowingMesh
                \brief Find all faces that share an edge with a given face
                \tparam Traits_ `CGAL::Region_growing::Region_growing_traits`
                \tparam Mesh_ Has model `CGAL::Surface_mesh`
                \cgalModels `RegionGrowingConnectivity`
            */

            template <class Traits_, class Mesh_>
            class Mesh_connectivity {

            public:
                using Face                    = typename Traits_::Element;
                ///< Must be equivalent to Mesh_::Face_index

                using Face_range              = CGAL::Iterator_range<CGAL::Face_around_face_iterator<Mesh_> >;
                ///< Result type of CGAL::faces_around_face()
                
                using Element_map             = typename Traits_::Element_map;
                ///< An `LvaluePropertyMap` that maps to a Face_index
                
                using Element_with_properties = typename Element_map::key_type;
                ///< Value type of the iterator in the input range

                /*!
                    The constructor requires a mesh to perform face searching process
                */
                Mesh_connectivity(const Mesh_& mesh) :
                m_mesh(mesh) { }

                /*!
                    From a query element `ewp`, this function retrieves the face via the element map and uses `CGAL::faces_around_face` to get a list of neighbor faces. The result is returned in `neighbors`.
                    \tparam Neighbors_ CGAL::Region_growing::Generalized_region_growing::Neighbors
                */
                template < class Neighbors_ >
                void get_neighbors(const Element_with_properties& ewp, Neighbors_& neighbors) {
                    Face face = get(m_elem_map, ewp);
                    neighbors.clear();
                    const Face_range& tmp = faces_around_face(m_mesh.halfedge(face), m_mesh);
                    for (typename Face_range::iterator it = tmp.begin(); it != tmp.end(); ++it) {
                        if ((*it) != Face(4294967295))
                            neighbors.push_back(*it);
                    }
                }

            private:
                const Mesh_& m_mesh;
                const Element_map m_elem_map = Element_map();
            };

        }
    }
}

#endif