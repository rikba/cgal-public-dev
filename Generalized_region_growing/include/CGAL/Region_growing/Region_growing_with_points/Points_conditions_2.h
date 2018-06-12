#ifndef GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H
#define GENERALIZED_REGION_GROWING_POINTS_CONDITIONS_2_H

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/squared_distance_2.h>

namespace CGAL {
    namespace Region_growing {
        namespace Region_growing_with_points {
            template<class Traits, class NormalMap>
            class Points_conditions_2 {
            public:
                using Kernel                  = typename Traits::Kernel;
                using Element_map             = typename Traits::Element_map;

                using Normal_map              = NormalMap;

                using Element_with_properties = typename Element_map::key_type;

                using Point_2  = typename Kernel::Point_2;
                using Line_2   = typename Kernel::Line_2;
                using Vector_2 = typename Kernel::Vector_2;
                using FT       = typename Kernel::FT;

                Points_conditions_2(const FT &epsilon, const FT &normal_threshold, const int &min_region_size) :
                        m_epsilon(epsilon),
                        m_normal_threshold(normal_threshold),
                        m_min_region_size(min_region_size) {}

                // Local condition
                template<class Region>
                bool is_in_same_region(const Element_with_properties &assigned_element,
                                       const Element_with_properties &unassigned_element,
                                       const Region &region) {

                    update(assigned_element, region);

                    Point_2 point_unassigned = get(m_elem_map, unassigned_element);

                    Vector_2 normal = get(m_normal_map, unassigned_element);
                    const FT normal_length = CGAL::sqrt(normal.squared_length());

                    Vector_2 normal_unassigned = normal / normal_length;

                    const FT distance_to_fit_line = CGAL::sqrt(
                            CGAL::squared_distance(point_unassigned, m_line_of_best_fit));
                    const FT cos_angle = CGAL::abs(normal_unassigned * m_normal_of_best_fit);

                    return (distance_to_fit_line <= m_epsilon && cos_angle >= m_normal_threshold);
                }

                // Global condition
                template<class Region>
                bool is_valid(const Region &region) const {
                    return (region.size() > m_min_region_size);
                }

                // Update the line of best fit
                template<class Region>
                void update(const Element_with_properties &assigned_element_with_properties, const Region &region) {

                    CGAL_precondition(region.size() > 0);

                    if (region.size() == 1) {
                        // The only point in the region is indeed `assigned_element_with_properties`
                        // The best fit line will be a line through this point with its normal being the point's normal

                        Point_2 point = get(m_elem_map, assigned_element_with_properties);
                        Vector_2 normal = get(m_normal_map, assigned_element_with_properties);
                        const FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_line_of_best_fit = Line_2(point, normal).perpendicular(point);
                        m_normal_of_best_fit = normal / normal_length;

                    } else {

                        // Extract the geometric Element (Point_2)
                        std::vector<Point_2> points;
                        for (typename Region::const_iterator it = region.begin(); it != region.end(); ++it)
                            points.push_back(get(m_elem_map, *it));

                        // Fit the region to a line
                        linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit,
                                                       CGAL::Dimension_tag<0>());

                        Vector_2 normal = m_line_of_best_fit.perpendicular(Point_2(0, 0)).to_vector();
                        const FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

            private:
                const Normal_map m_normal_map = Normal_map();
                const Element_map m_elem_map = Element_map();
                Line_2 m_line_of_best_fit;
                Vector_2 m_normal_of_best_fit;
                const FT &m_epsilon;
                const FT &m_normal_threshold;
                const int &m_min_region_size;
            };
        }
    }
}
#endif
