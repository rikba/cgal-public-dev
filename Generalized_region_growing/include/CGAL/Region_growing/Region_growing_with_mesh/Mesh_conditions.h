#ifndef CGAL_GRG_MESH_CONDITIONS_H
#define CGAL_GRG_MESH_CONDITIONS_H

#include <Surface_mesh.h>
#include <Iterator_range.h>
#include <CGAL/Cartesian_converter.h>


namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            template <class Traits, class Mesh>
            class Mesh_conditions {

                template<class ...>
                using void_t = void;

                template<class Kernel, class = void>
                class Get_sqrt {

                    typedef typename Kernel::FT FT;

                public:
                    FT operator()(const FT &value) const {

                        return FT(CGAL::sqrt(CGAL::to_double(value)));
                    }
                };

                template<class Kernel>
                class Get_sqrt<Kernel, void_t<typename Kernel::Sqrt> > : Kernel::Sqrt { };

            public:
                using Kernel = typename Traits::Kernel;
                using Face_index = typename Mesh::Face_index;
                using Vertex_index = typename Mesh::Vertex_index;
                using Vertex_range = typename Mesh::Vertex_range;
                using Plane_3 = typename Kernel::Plane_3;
                using Point_3 = typename Kernel::Point_3;
                using Element_map = typename Traits::Element_map;
                using FT = typename Kernel::FT;
                using Sqrt                    = Get_sqrt<Kernel>;

                Meshes_conditions(const Mesh& mesh, const FT& epsilon, const FT& normal_threshold) :
                m_mesh(mesh),
                m_epsilon(epsilon),
                m_normal_threshold(normal_threshold) { }
                
                template <class Region, class ElementWithProperties>
                bool is_in_same_region(const ElementWithProperties& assigned_element, const ElementWithProperties& unassigned_element, const Region& region) {

                    Face_index face = get(m_elem_map, unassigned_element);

                    const Point_3 face_centroid = get_centroid(face);
                    const Vector_3 face_normal = get_normal(face);

                    const Local_FT distance_to_fit_plane = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(face_centroid), m_plane_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(face_normal) * m_normal_of_best_fit);

                    return (distance_to_fit_plane <= m_epsilon && cos_angle >= m_normal_threshold);

                }

                template <class Region>
                bool update_shape(const Region& region) {

                    CGAL_precondition(region.end() - region.begin() != 0);

                    if (region.end() - region.begin() == 1) {
                        
                        Face_index face = get(m_elem_map, *region.begin());

                        const Point_3 face_centroid = get_centroid(face);
                        const Vector_3 face_normal = get_normal(face);

                        const FT normal_length = m_sqrt(face_normal.squared_length());

                        m_plane_of_best_fit = m_to_local_converter(Plane_3(face_centroid, face_normal));
                        m_normal_of_best_fit = m_to_local_converter(face_normal / normal_length);

                    } else {

                        // Extract the points of the region (Point_3)
                        std::vector<Local_point_3> points;
                        // The region is formed by face indices
                        for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it) {
                            
                            // Get vertices of each face and push them to `points`
                            Vertex_range vr = vertices_around_face(m_mesh.halfedge(*it), m_mesh);
                            for (typename Vertex_range::const_iterator vit = vr.begin(); vit != vr.end(); ++vit)
                                points.push_back(m_to_local_converter(get(m_elem_map, *it)));

                        }

                        // Fit the region to a plane
                        linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, CGAL::Dimension_tag<0>());

                        Local_vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
                        const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

            private:

                Point_3 get_centroid(const Face_index& face) {
                    Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);

                    // Calculate centroid
                    Point_3 centroid = Point_3(0, 0, 0);
                    int n = 0;
                    for (Vertex_range::iterator it = vr.begin(); it != vr.end(); ++it, ++n)
                        centroid += point(*it);
                    CGAL_precondition(n > 2);
                    centroid /= n;

                    return centroid;
                }

                Vector_3 get_normal(const Face_index& face) {
                    Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);
                    return CGAL::normal(point(*vr.begin()), point(*(vr.begin()+1)), point(*(vr.begin()+2)));
                }

                const Mesh& m_mesh;
                const Element_map m_elem_map = Element_map();
                const FT& m_epsilon;
                const FT& m_normal_threshold;
                const Sqrt                    m_sqrt;
                const To_local_converter      m_to_local_converter;
                const To_input_converter      m_to_input_converter;
                Local_plane_3                 m_plane_of_best_fit;
                Local_vector_3                m_normal_of_best_fit;
            };

        }
    }
}

#endif