#ifndef CGAL_REGION_GROWING_POINTS_CONDITIONS_2_H
#define CGAL_REGION_GROWING_POINTS_CONDITIONS_2_H

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {

    namespace Region_growing {

        namespace Region_growing_with_points {

            /*! 
                \ingroup PkgGeneralizedRegionGrowingPoints
                \brief Local and global conditions for the region growing algorithm on 3D point cloud.
                \tparam Traits_ CGAL::Region_growing::Region_growing_traits.
                \tparam NormalMap An `LvaluePropertyMap` that maps to a vector, representing the normal associated with the point.
                \cgalModels `RegionGrowingConditions`
            */

            template<class Traits_, class NormalMap>
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
                using Kernel                  = typename Traits_::Kernel;
                using Element_map             = typename Traits_::Element_map;
                using Normal_map              = NormalMap;

                using Element_with_properties = typename Element_map::key_type;
                ///< Value type of the iterator in the input range

                using Point_2                 = typename Kernel::Point_2; ///< Point type
                using Line_2                  = typename Kernel::Line_2; ///< Line type
                using Vector_2                = typename Kernel::Vector_2; ///< Vector type
                using FT                      = typename Kernel::FT; ///< Number type

                ///< \cond SKIP_IN_MANUAL
                using Sqrt                    = Get_sqrt<Kernel>;
                using Local_kernel            = Exact_predicates_inexact_constructions_kernel;
                using To_local_converter      = Cartesian_converter<Kernel, Local_kernel>;
                using Local_point_2           = Local_kernel::Point_2;
                using Local_line_2            = Local_kernel::Line_2;
                using Local_vector_2          = Local_kernel::Vector_2;
                using Local_FT                = Local_kernel::FT;
                ///< \endcond

                /*!
                    Each region is represented by a line. The constructor requires three parameters, in order: the maximum distance from a point to the region, the minimum dot product between the normal associated with the point and the normal of the region, and the minimum number of points a region must have.
                */
                Points_conditions_2(const FT epsilon, const FT normal_threshold, const size_t min_region_size) :
                m_epsilon(epsilon),
                m_normal_threshold(normal_threshold),
                m_min_region_size(min_region_size),
                m_sqrt(Sqrt()) {}

                /*!
                    Local condition that checks if a new point in `unassigned_element` is similar to the point `assigned_element` and its enclosing region. Each item in `Region_` has the same structure as in the input range, i.e. has the type `Element_with_properties`.
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */
                template < class Region_ >
                bool is_in_same_region(const Element_with_properties &assigned_element,
                                       const Element_with_properties &unassigned_element,
                                       const Region_ &region) {

                    const Point_2& point_unassigned = get(m_elem_map, unassigned_element);
                    const Vector_2& normal = get(m_normal_map, unassigned_element);

                    const FT normal_length = m_sqrt(normal.squared_length());
                    Vector_2 normal_unassigned = normal / normal_length;

                    // Must use Local_FT because fit line is of local kernel
                    const Local_FT distance_to_fit_line = CGAL::sqrt(CGAL::squared_distance(m_to_local_converter(point_unassigned), m_line_of_best_fit));
                    const Local_FT cos_angle = CGAL::abs(m_to_local_converter(normal_unassigned) * m_normal_of_best_fit);

                    return (distance_to_fit_line <= m_epsilon && cos_angle >= m_normal_threshold);
                }

                /*!
                    Global condition that checks if a region size is large enough to be accepted
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */
                template < class Region_ >
                inline bool is_valid(const Region_ &region) const {
                    return (region.size() >= m_min_region_size);
                }

                /*!
                    Update the class' best fit plane that will be used later by the local condition.
                    \tparam Region_ CGAL::Region_growing::Generalized_region_growing::Region
                */
                template < class Region_ >
                void update_shape(const Region_ &region) {

                    CGAL_precondition(region.end() != region.begin());

                    if (region.size() == 1) {
                        // The best fit line will be a line through this point with its normal being the point's normal

                        const Point_2& point = get(m_elem_map, *region.begin());
                        const Vector_2& normal = get(m_normal_map, *region.begin());
                        const FT normal_length = m_sqrt(normal.squared_length());

                        m_line_of_best_fit = m_to_local_converter(Line_2(point, normal).perpendicular(point));
                        m_normal_of_best_fit = m_to_local_converter(normal / normal_length);

                    } else {

                        // Extract the geometric Element (Point_2)
                        int i = 0;
                        std::vector<Local_point_2> points(region.size());
                        for (typename Region_::const_iterator it = region.begin(); it != region.end(); ++it, ++i)
                            points[i] = m_to_local_converter(get(m_elem_map, *it));

                        // Fit the region to a line
                        Local_point_2 centroid; // unused

                        #ifndef CGAL_EIGEN2_ENABLED
                            linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Kernel(), Default_diagonalize_traits<typename Kernel::FT, 2>());
                        #else 
                            linear_least_squares_fitting_2(points.begin(), points.end(), m_line_of_best_fit, centroid, CGAL::Dimension_tag<0>(), Kernel(), Eigen_diagonalize_traits<typename Kernel::FT, 2>());
                        #endif

                        Local_vector_2 normal = m_line_of_best_fit.perpendicular(m_line_of_best_fit.point(0)).to_vector();
                        const Local_FT normal_length = CGAL::sqrt(normal.squared_length());

                        m_normal_of_best_fit = normal / normal_length;
                    }
                }

            private:
                const Normal_map                m_normal_map = Normal_map();
                const Element_map               m_elem_map = Element_map();
                const FT                        m_epsilon;
                const FT                        m_normal_threshold;
                const size_t                    m_min_region_size;
                const Sqrt                      m_sqrt;
                const To_local_converter        m_to_local_converter;
                Local_line_2                    m_line_of_best_fit;
                Local_vector_2                  m_normal_of_best_fit;
            };

        } // namespace Region_growing_with_points

    } // namespace Region_growing

} // namespace CGAL
#endif
