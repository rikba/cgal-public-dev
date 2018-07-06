#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            template<class Traits, class NormalMap>
            class Points_conditions_2 {

                template<class ...>
                using void_t = void;

                template<class Kernel, class = void>
                class Get_sqrt {
                    typedef typename Kernel::FT FT;
                public:
                    FT operator()(const FT &value) const {
                        return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
                    }
                };

                template<class Kernel>
                class Get_sqrt<Kernel, void_t<typename Kernel::Sqrt> > : Kernel::Sqrt { };

            public:
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;
                using Normal_map              = NormalMap;

                using Element_with_properties = typename Element_map::key_type;
                using Point_2                 = typename Kernel::Point_2;
                using Line_2                  = typename Kernel::Line_2;
                using Vector_2                = typename Kernel::Vector_2;
                using FT                      = typename Kernel::FT;

                using Sqrt                    = Get_sqrt<Kernel>;

                using Local_kernel            = Exact_predicates_inexact_constructions_kernel;
                using To_local_converter      = Cartesian_converter<Kernel, Local_kernel>;
                using Local_point_2           = Local_kernel::Point_2;
                using Local_line_2            = Local_kernel::Line_2;
                using Local_vector_2          = Local_kernel::Vector_2;
                using Local_FT                = Local_kernel::FT;


                Points_conditions_2(const FT &epsilon, const FT &normal_threshold, const int &min_region_size) :
                    m_epsilon(epsilon),
                    m_normal_threshold(normal_threshold),
                    m_min_region_size(min_region_size),
                    m_sqrt(Sqrt()) {}

                // Local condition
                template<class Region>
                bool is_in_same_region(const Element_with_properties &assigned_element,
                                       const Element_with_properties &unassigned_element,
                                       const Region &region) {

                    Point_2 point_unassigned = get(m_elem_map, unassigned_element);
                    Vector_2 normal = get(m_normal_map, unassigned_element);

                    const FT normal_length = m_sqrt(normal.squared_length());
                    Vector_2 normal_unassigned = normal / normal_length;

                    // Must use Local_FT because fit line is of local kernel
                    const Local_FT distance_to_fit_line = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(point_unassigned), m_line_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(normal_unassigned) * m_normal_of_best_fit);

                    return (distance_to_fit_line <= m_epsilon && cos_angle >= m_normal_threshold);
                }

                // Global condition
                template<class Region>
                bool is_valid(const Region &region) const {

                    return (region.size() >= m_min_region_size);
                }

                // Update the line of best fit
                template<class Region>
                void update_shape(const Region &region) {

                    CGAL_precondition(region.end() - region.begin() != 0);

                    if (region.end() - region.begin() == 1) {
                        // The best fit line will be a line through this point with its normal being the point's normal

                        Point_2 point = get(m_elem_map, *region.begin());
                        Vector_2 normal = get(m_normal_map, *region.begin());
                        const FT normal_length = m_sqrt(normal.squared_length());

                        m_line_of_best_fit = m_to_local_converter(Line_2(point, normal).perpendicular(point));
                        m_normal_of_best_fit = m_to_local_converter(normal / normal_length);

                    } else {

                        // Extract the geometric Element (Point_2)
                        std::vector<Local_point_2> points;
                        for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it)
                            points.push_back(m_to_local_converter(get(m_elem_map, *it)));

                        // Fit the region to a line
                        linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit, CGAL::Dimension_tag<0>());

                        Local_vector_2 normal = m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
                        const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

            private:
                const Normal_map m_normal_map = Normal_map();
                const Element_map m_elem_map = Element_map();
                const FT &m_epsilon;
                const FT &m_normal_threshold;
                const int &m_min_region_size;
                const Sqrt m_sqrt;
                const To_local_converter m_to_local_converter;
                Local_line_2 m_line_of_best_fit;
                Local_vector_2 m_normal_of_best_fit;
            };

        } // namespace Region_growing_with_points

    } // namespace Region_growing

} // namespace CGAL
#endif
