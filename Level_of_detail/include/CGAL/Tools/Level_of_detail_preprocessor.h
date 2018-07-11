#ifndef CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <vector>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/constructions_d.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/property_map.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include "../../../test/Level_of_detail/include/Mylog.h"
#include <CGAL/Tools/Level_of_detail_projector.h>
#include <CGAL/Tools/Level_of_detail_utils.h>
#include <CGAL/Tools/Level_of_detail_interior_boundary_extractor.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_preprocessor {

		public:
			typedef KernelTraits   Traits;
			typedef InputContainer Container;

			typedef typename Traits::FT 	 FT;
			typedef typename Traits::Point_2 Point_2;
			typedef typename Traits::Point_3 Point_3;
			typedef typename Traits::Plane_3 Plane_3;

			typedef CGAL::Search_traits_3<Traits>                     Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                 Fuzzy_sphere;
			typedef typename Neighbor_search::Tree 					  Tree;

			typedef std::map<int, Point_2> 												  Projected_points;
			typedef Level_of_detail_simple_projector<Traits, Container, Projected_points> Ground_projector;

			using Index          = int;
			using Const_iterator = typename Container::const_iterator;
			using Index_map      = typename Container:: template Property_map<Index>;
			using Log 			 = CGAL::LOD::Mylog;

			Level_of_detail_preprocessor() : m_scale(FT(1)), m_alpha(-FT(1)), m_use_alpha_shapes(false), m_silent(false) { }

			void make_silent(const bool new_state) {
				m_silent = new_state;
			}

			void set_scale(const FT new_value) {
				
				assert(new_value >= FT(0));
				m_scale = new_value;
			}

			void set_alpha(const FT new_value) {

				assert(new_value > FT(0));
				m_alpha = new_value;
			}

			void use_alpha_shapes(const bool new_state) {
				m_use_alpha_shapes = new_state;
			}

			template<class Planes>
			int get_planes(const Container &input, Planes &planes) {

				auto number_of_planes = -1;
				create_indices(input);
				
				planes.clear();

				for (Const_iterator it = input.begin(); it != input.end(); ++it)
					if (m_indices[*it] >= 0) 
						planes[m_indices[*it]].push_back(*it);

				number_of_planes = planes.size();
				return number_of_planes;
			}

			// Later I can adapt boundary_clutter to some other data structure where I also take
			// into account the type of the clutter: undetected, cylinder, sphere and so on. The type is saved in property_map<Types>.
			template<class Indices, class Boundary_data>
			int get_boundary_points(
				const Container &input, const Indices &boundary_mapping, const Indices &interior_mapping, const bool with_shape_detection, 
				Boundary_data &building_boundaries, Boundary_data &boundary_clutter) {

				boundary_clutter.clear();
				create_indices(input);

				if (m_use_alpha_shapes) add_interior_boundary_to_clutter(input, boundary_mapping, boundary_clutter);
				else add_boundary_points(boundary_mapping, with_shape_detection, building_boundaries, boundary_clutter);

				add_interior_boundary_to_clutter(input, interior_mapping, boundary_clutter);
				return 0;
			}

			// This is very slow algorithm. Should be improved later.
			template<class Projected_points, class Point_sets>
			int clean_projected_points(Projected_points &projected_points, Point_sets &point_sets) {

				using Traits1   = CGAL::Simple_cartesian<double>;
				using Point_3ft = Traits1::Point_3;

				if (projected_points.size() == 0) return 0;
				assert(!point_sets.empty());

				Projected_points cleaned_points;
				Point_sets cleaned_sets;

				const auto num_neighbours = 2;

				CGAL::Identity_property_map<Point_3ft> pmap1;
				for (typename Point_sets::const_iterator it = point_sets.begin(); it != point_sets.end(); ++it) {
					const auto set_index = (*it).first;

					const size_t num_points = (*it).second.size();
					std::vector<Point_3ft> points1(num_points);
					std::vector<Point_3> points2(num_points);

					for (size_t i = 0; i < num_points; ++i) {
						const auto point_index = (*it).second[i];

						const Point_2 &p = projected_points.at(point_index);

						const double x = CGAL::to_double(p.x());
						const double y = CGAL::to_double(p.y());
						const double z = 0.0;

						points1[i] = Point_3ft(x, y, z);
						points2[i] = Point_3(p.x(), p.y(), FT(0));
					}
						
					const FT average_spacing = static_cast<FT>(CGAL::compute_average_spacing<CGAL::Sequential_tag>(points1.begin(), points1.end(), pmap1, num_neighbours, Traits1()));
					Tree tree(points2.begin(), points2.end());

					std::vector<Point_3> qs;
					for (size_t i = 0; i < num_points; ++i) {
						const auto point_index = (*it).second[i];
						
						const Point_3 &p = points2[i];
						const Fuzzy_sphere sphere(p, m_scale * average_spacing);

						qs.clear();	
						tree.search(std::back_inserter(qs), sphere);

						if (qs.empty() || (qs.size() == 1 && qs[0] == p)) continue;

						cleaned_points[point_index] = Point_2(p.x(), p.y());
						cleaned_sets[set_index].push_back(point_index);
					}
				}

				const auto number_of_outliers = static_cast<int>(projected_points.size() - cleaned_points.size());

				projected_points = cleaned_points;
				point_sets = cleaned_sets;

				return number_of_outliers;
			}

		private:
			Index_map m_indices;

			FT m_scale;
			FT m_alpha;

			bool m_use_alpha_shapes;
			bool m_silent;

			void create_indices(const Container &input) {
				boost::tie(m_indices, boost::tuples::ignore) = input. template property_map<Index>("index");
			}

			template<class Indices, class Boundary_data>
			int add_boundary_points(const Indices &mapping, const bool with_shape_detection, Boundary_data &building_boundaries, Boundary_data &boundary_clutter) {

				auto number_of_boundaries = -1;
				building_boundaries.clear();

				if (mapping.empty()) return number_of_boundaries;

				for (size_t i = 0; i < mapping.size(); ++i) {
					if (m_indices[mapping[i]] >= 0 && with_shape_detection) building_boundaries[m_indices[mapping[i]]].push_back(mapping[i]);
					else boundary_clutter[0].push_back(mapping[i]);
				}

				number_of_boundaries = building_boundaries.size();
				return number_of_boundaries;
			}


			template<class Indices, class Boundary_data>
			void add_interior_boundary_to_clutter(const Container &input, const Indices &interior_mapping, Boundary_data &boundary_clutter) {
				if (interior_mapping.empty()) return;

				// (1) Project points onto the ground.
				Projected_points building_interior_projected;
				project_points_onto_ground<Indices, Boundary_data>(input, interior_mapping, building_interior_projected);


				// (2) Extract boundaries.
				Level_of_detail_interior_boundary_extractor<Traits, Container, Indices> extractor;

				Indices result;
				extractor.set_alpha(m_alpha);
				extractor.make_silent(m_silent);
				extractor.extract(input, interior_mapping, building_interior_projected, result);


				// (3) Add extracted points to the clutter.
				Boundary_data stub;
				add_boundary_points(result, false, stub, boundary_clutter);
			}

			template<class Indices, class Boundary_data>
			void project_points_onto_ground(const Container &input, const Indices &interior_mapping, Projected_points &building_interior_projected) {
				if (interior_mapping.empty()) return;

				Boundary_data building_interior, stub;
				add_boundary_points(interior_mapping, false, stub, building_interior);

				Plane_3 base_ground_plane(FT(0), FT(0), FT(1), FT(0));
				building_interior_projected.clear();
				
				Ground_projector projector; 
				projector.project(input, building_interior, base_ground_plane, building_interior_projected);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H