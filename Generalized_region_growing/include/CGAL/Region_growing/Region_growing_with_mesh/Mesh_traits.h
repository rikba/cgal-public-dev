#ifndef CGAL_GRG_MESH_TRAITS_H
#define CGAL_GRG_MESH_TRAITS_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            template<class InputRange, class ElementMap, class K>
            class Mesh_traits {
            public:
                using Kernel                  = K;
                using Input_range             = InputRange;
                using Element_map             = ElementMap;
                using Element                 = typename Element_map::value_type;
                // Element is a face of a mesh, which can be CGAL::Surface_mesh<K::Point_3>::Face_index
            };

        }
    }
}

#endif