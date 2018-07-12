#ifndef CGAL_GRG_POINTS_CONDITIONS_2_H
#define CGAL_GRG_POINTS_CONDITIONS_2_H

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
                using Input_range             = typename Traits::Input_range;
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


                Points_conditions_2(const Input_range& input_range, const FT &epsilon, const FT &normal_threshold, const int &min_region_size) :
                    m_input_range(input_range),
                    m_epsilon(epsilon),
                    m_normal_threshold(normal_threshold),
                    m_min_region_size(min_region_size),
                    m_sqrt_object(Sqrt()) {}

                // Local condition
                template < class Region_ >
                bool is_in_same_region(const int &assigned_element, const int &unassigned_element,const Region_ &region) {

                    Point_2 point_unassigned = get(m_elem_map, get_data_from_index(unassigned_element));
                    Vector_2 normal = get(m_normal_map, get_data_from_index(unassigned_element));

                    const FT normal_length = m_sqrt_object(normal.squared_length());
                    Vector_2 normal_unassigned = normal / normal_length;

                    // Must use Local_FT because fit line is of local kernel
                    const Local_FT distance_to_fit_line = m_sqrt_object(CGAL::squared_distance(m_to_local_converter(point_unassigned), m_line_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(normal_unassigned) * m_normal_of_best_fit);

                    return (distance_to_fit_line <= m_epsilon && cos_angle >= m_normal_threshold);
                }

                // Global condition
                template < class Region_ >
                inline bool is_valid(const Region_ &region) const {
                    return (region.end() - region.begin() >= m_min_region_size);
                }

                // Update the line of best fit
                template < class Region_ >
                void update_shape(const Region_ &region) {

                    CGAL_precondition(region.end() != region.begin());

                    if (region.end() - region.begin() == 1) {
                        // The best fit line will be a line through this point with its normal being the point's normal

                        Point_2 point = get(m_elem_map, get_data_from_index(*region.begin()));
                        Vector_2 normal = get(m_normal_map, get_data_from_index(*region.begin()));
                        const FT normal_length = m_sqrt_object(normal.squared_length());

                        m_line_of_best_fit = m_to_local_converter(Line_2(point, normal).perpendicular(point));
                        m_normal_of_best_fit = m_to_local_converter(normal / normal_length);

                    } else {

                        // Extract the geometric Element (Point_2)
                        int i = 0;
                        std::vector<Point_2> points(region.end()-region.begin());
                        for (typename Region_::const_iterator it = region.begin(); it != region.end(); ++it, ++i)
                            points[i] = (m_to_local_converter(get(m_elem_map, get_data_from_index(*it))));

                        // Fit the region to a line
                        Line_2 temp_best_fit;
                        linear_least_squares_fitting_2(points.begin(), points.end(), temp_best_fit, CGAL::Dimension_tag<0>());

                        m_line_of_best_fit = m_to_local_converter(temp_best_fit);

                        Local_vector_2 normal = m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
                        const Local_FT normal_length = m_sqrt_object(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

                inline Element_with_properties get_data_from_index(const int& i) const {
                    return *(m_input_range.begin() + i);
                }

            private:
                const Input_range &             m_input_range;
                const Normal_map                m_normal_map = Normal_map();
                const Element_map               m_elem_map = Element_map();
                const FT &                      m_epsilon;
                const FT &                      m_normal_threshold;
                const int &                     m_min_region_size;
                const Sqrt                      m_sqrt_object;
                const To_local_converter        m_to_local_converter;
                Local_line_2                    m_line_of_best_fit;
                Local_vector_2                  m_normal_of_best_fit;
            };

        } // namespace Region_growing_with_points

    } // namespace Region_growing

} // namespace CGAL
#endif
