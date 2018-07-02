#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_3_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_3_H

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            template<class Traits, class NormalMap>
            class Points_conditions_3 {

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
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;
                using Normal_map              = NormalMap;

                using Element_with_properties = typename Element_map::key_type;
                using Point_3                 = typename Kernel::Point_3;
                using Plane_3                 = typename Kernel::Plane_3;
                using Vector_3                = typename Kernel::Vector_3;
                using FT                      = typename Kernel::FT;

                using Sqrt                    = Get_sqrt<Kernel>;

                using Local_kernel            = Exact_predicates_inexact_constructions_kernel;
                using To_local_converter      = Cartesian_converter<Kernel, Local_kernel>;
                using To_input_converter      = Cartesian_converter<Local_kernel, Kernel>;
                using Local_point_3           = Local_kernel::Point_3;
                using Local_plane_3           = Local_kernel::Plane_3;
                using Local_vector_3          = Local_kernel::Vector_3;
                using Local_FT                = Local_kernel::FT;

                Points_conditions_3(const FT &epsilon, const FT &normal_threshold, const int &min_region_size) :
                    m_epsilon(epsilon),
                    m_normal_threshold(normal_threshold),
                    m_min_region_size(min_region_size),
                    m_sqrt(Sqrt()) {}

                // Local condition
                template<class Region>
                bool is_in_same_region(const Element_with_properties &assigned_element,
                                       const Element_with_properties &unassigned_element,
                                       const Region &region) {

                    Point_3 point_unassigned = get(m_elem_map, unassigned_element);
                    Vector_3 normal = get(m_normal_map, unassigned_element);

                    const FT normal_length = m_sqrt(normal.squared_length());
                    Vector_3 normal_unassigned = normal / normal_length;

                    // Must use Local_FT because the fit plane is of local kernel
                    const Local_FT distance_to_fit_plane = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(point_unassigned), m_plane_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(normal_unassigned) * m_normal_of_best_fit);

                    return (distance_to_fit_plane <= m_epsilon && cos_angle >= m_normal_threshold);
                }

                // Global condition
                template<class Region>
                bool is_valid(const Region &region) const {

                    return (region.size() > m_min_region_size);
                }

                // Update the plane of best fit
                template<class Region>
                void update_shape(const Region &region) {

                    CGAL_precondition(region.end() - region.begin() != 0);

                    if (region.end() - region.begin() == 1) {
                        // The only point in the region is indeed `assigned_element_with_properties`
                        // The best fit plane will be a plane through this point with its normal being the point's normal

                        Point_3 point = get(m_elem_map, *region.begin());
                        Vector_3 normal = get(m_normal_map, *region.begin());
                        const FT normal_length = m_sqrt(normal.squared_length());

                        m_plane_of_best_fit = m_to_local_converter(Plane_3(point, normal));
                        m_normal_of_best_fit = m_to_local_converter(normal / normal_length);

                    } else {

                        // Extract the geometric Element (Point_3)
                        std::vector<Local_point_3> points;
                        for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it)
                            points.push_back(m_to_local_converter(get(m_elem_map, *it)));

                        // Fit the region to a plane
                        linear_least_squares_fitting_3(points.begin(), points.end(), m_plane_of_best_fit, CGAL::Dimension_tag<0>());

                        Local_vector_3 normal = m_plane_of_best_fit.orthogonal_vector();
                        const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

                Plane_3 plane_of_best_fit() const {

                    return m_to_input_converter(m_plane_of_best_fit);
                }

                Vector_3 normal_of_best_fit() const {

                    return m_to_input_converter(m_normal_of_best_fit);
                }

            private:
                const Normal_map              m_normal_map = Normal_map();
                const Element_map             m_elem_map = Element_map();
                const FT &                    m_epsilon;
                const FT &                    m_normal_threshold;
                const int &                   m_min_region_size;
                const Sqrt                    m_sqrt;
                const To_local_converter      m_to_local_converter;
                const To_input_converter      m_to_input_converter;
                Local_plane_3                 m_plane_of_best_fit;
                Local_vector_3                m_normal_of_best_fit;
                size_t                        m_region_size = 0;
            }; // Points_conditions_3<Traits, NormalMap>

        }
    }
}
#endif
