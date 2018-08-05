#ifndef CGAL_REGION_GROWING_MESH_CONDITIONS_H
#define CGAL_REGION_GROWING_MESH_CONDITIONS_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Cartesian_converter.h>


namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_mesh {

            /*!
                \ingroup PkgGeneralizedRegionGrowingMesh
                \brief Local and global conditions for the region growing algorithm on a mesh.
                \tparam Traits_ `CGAL::Region_growing::Region_growing_traits`
                \tparam Mesh_ Has model `CGAL::Surface_mesh`
                \cgalModels `RegionGrowingConditions`
            */

            template <class Traits_, class Mesh_>
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
                using Kernel                  = typename Traits_::Kernel;
                using Face                    = typename Traits_::Element;
                ///< Must be equivalent to Mesh_::Face_index
                using Vertex_index            = typename Mesh_::Vertex_index;

                using Vertex_range            = Iterator_range<Vertex_around_face_iterator<Mesh_> >;
                ///< Result type of CGAL::vertices_around_face()

                using Plane_3                 = typename Kernel::Plane_3; ///< Plane type
                using Point_3                 = typename Kernel::Point_3; ///< Point type
                using Vector_3                = typename Kernel::Vector_3; ///< Vector type

                using Element_map             = typename Traits_::Element_map; 
                ///< An `LvaluePropertyMap` that maps to a Face_index

                using FT                      = typename Kernel::FT; ///< Number type

                #ifndef DOXYGEN_RUNNING
                    using Sqrt                    = Get_sqrt<Kernel>;
                    using Local_kernel            = Exact_predicates_inexact_constructions_kernel;
                    using To_local_converter      = Cartesian_converter<Kernel, Local_kernel>;
                    using Local_point_3           = Local_kernel::Point_3;
                    using Local_plane_3           = Local_kernel::Plane_3;
                    using Local_vector_3          = Local_kernel::Vector_3;
                    using Local_FT                = Local_kernel::FT;
                #endif

                using Element_with_properties = typename Element_map::key_type;
                ///< Value type of the iterator in the input range

                /*!
                    Each region is represented by a plane. The constructor requires three parameters, in order: the mesh, the maximum distance from a point to the region, and the minimum dot product between the normal associated with the point and the normal of the region.
                */                
                Mesh_conditions(const Mesh_& mesh, const FT epsilon, const FT normal_threshold) :
                m_mesh(mesh),
                m_epsilon(epsilon),
                m_normal_threshold(normal_threshold) { }
                
                /*!
                    Local condition that checks if a new face in `unassigned_element` is similar to the face `assigned_element` and its enclosing region `region`.
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */
                template < class Region_ >
                bool is_in_same_region(const Element_with_properties &assigned_element,
                                       const Element_with_properties &unassigned_element,
                                       const Region_ &region) {

                    const Face& face = get(m_elem_map, unassigned_element);

                    Point_3 face_centroid;
                    get_centroid(face, face_centroid);
                    Vector_3 face_normal;
                    get_normal(face, face_normal);

                    const Local_FT distance_to_fit_plane = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(face_centroid), m_plane_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(face_normal) * m_normal_of_best_fit);

                    return (distance_to_fit_plane <= m_epsilon && cos_angle >= m_normal_threshold);

                }

                /*!
                    Update the class' best fit plane that will be used later by the local condition.
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */

                template < class Region_ >
                void update_shape(const Region_& region) {

                    CGAL_precondition(region.size() != 0);

                    if (region.size() == 1) {
                        
                        const Face& face = get(m_elem_map, *region.begin());

                        Point_3 face_centroid;
                        get_centroid(face, face_centroid);
                        Vector_3 face_normal;
                        get_normal(face, face_normal);

                        const FT normal_length = m_sqrt(face_normal.squared_length());

                        m_plane_of_best_fit = m_to_local_converter(Plane_3(face_centroid, face_normal));
                        m_normal_of_best_fit = m_to_local_converter(face_normal / normal_length);

                    } else {

                        // Extract the points of the region (Point_3)
                        std::vector<Local_point_3> points;
                        // The region is formed by face indices
                        for (typename Region_::const_iterator it = region.begin(); it != region.end(); ++it) {

                            const Face& face = get(m_elem_map, *it);
                            
                            // Get vertices of each face and push them to `points`
                            Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);
                            for (typename Vertex_range::const_iterator vit = vr.begin(); vit != vr.end(); ++vit)
                                points.push_back(m_to_local_converter(m_mesh.point(*vit)));

                        }
                        
                        Point_3 centroid; // unused
                        // Fit the region to a plane
                        #ifndef CGAL_EIGEN3_ENABLED
                            linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Kernel(), Default_diagonalize_traits<typename Kernel::FT, 3>());
                        #else 
                            linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Kernel(), Eigen_diagonalize_traits<typename Kernel::FT, 3>());
                        #endif

                        Local_vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
                        const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

                /*!
                    Global condition that checks the validity of a region. Always return true if the region is not empty.
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */
                template < class Region_ >
                inline bool is_valid(const Region_& region) {
                    return region.size() > 0;
                }

            private:

                void get_centroid(const Face& face, Point_3& centroid) {
                    Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);

                    // Calculate centroid
                    Vector_3 centr = Vector_3(0, 0, 0); // temporarily set as a vector to perform + and / operators
                    int n = 0;
                    for (typename Vertex_range::iterator it = vr.begin(); it != vr.end(); ++it, ++n) {
                        Point_3 tmp = m_mesh.point(*it);
                        centr += Vector_3(tmp.x(), tmp.y(), tmp.z());
                    }
                    CGAL_precondition(n > 2);
                    centr /= n;

                    centroid = Point_3(centr.x(), centr.y(), centr.z()); // convert back to point
                }

                void get_normal(const Face& face, Vector_3& normal) {
                    Vertex_range vr = vertices_around_face(m_mesh.halfedge(face), m_mesh);
                    typename Vertex_range::iterator it = vr.begin();
                    normal = CGAL::normal(m_mesh.point(*it++), m_mesh.point(*it++), m_mesh.point(*it++));
                }

                const Mesh_&                  m_mesh;
                const Element_map             m_elem_map = Element_map();
                const FT                      m_epsilon;
                const FT                      m_normal_threshold;
                const Sqrt                    m_sqrt;
                const To_local_converter      m_to_local_converter;
                Local_plane_3                 m_plane_of_best_fit;
                Local_vector_3                m_normal_of_best_fit;
            };

        }
    }
}

#endif