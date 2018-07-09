#ifndef GENERALIZED_REGION_GROWING_MESHES_TRAITS_H
#define GENERALIZED_REGION_GROWING_MESHES_TRAITS_H

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_meshes {

            template<class InputRange, class ElementMap, class K>
            class Meshes_traits {
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