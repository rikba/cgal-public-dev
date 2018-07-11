#ifndef CGAL_LEVEL_OF_DETAIL_GRID_SIMPLIFY_H
#define CGAL_LEVEL_OF_DETAIL_GRID_SIMPLIFY_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#define PN "\r\n"
#else 
#define PSR "/" 
#define PN "\n"
#endif

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <cmath>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include "../../../../../../test/Level_of_detail/include/debugging/Mylog.h"
#include "../../../../../../test/Level_of_detail/include/temporary/Level_of_detail_enumerations.h"

namespace CGAL {

	namespace Level_of_detail {

		// WARNING: DOES NOT WORK IF LATER I NEED TO REUSE ANY DATA FROM THE ORIGINAL INPUT LIKE NORMALS E.G.
		// ONLY THE CLOSEST NEW POINT METHOD WILL WORK IF REFACTORED PROPERLY.
		// IN ALL OTHER CASES I INTRODUCE NEW POINTS, WHICH DO NOT EXIST IN THE ORIGINAL INPUT!

		// Grid based simplification.
		template<class KernelTraits, class BoundaryData, class ProjectedPoints>
		class Level_of_detail_grid_simplify {

		public:
			// Fields.
			typedef KernelTraits 	Kernel;
			typedef BoundaryData 	Boundary_data;
			typedef ProjectedPoints Projected_points;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_2 Point_2;
			typedef typename Kernel::Line_2  Line_2;

			typedef std::pair<int, Point_2> 					 				Projected_point;
			typedef typename CGAL::Second_of_pair_property_map<Projected_point> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					 Search_traits_2;
			typedef CGAL::Search_traits_adapter<Projected_point, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 					     Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					 Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                         Fuzzy_tree;

			typedef typename Boundary_data::const_iterator    Boundary_iterator;
			typedef typename Projected_points::const_iterator Point_iterator;
			typedef typename Neighbor_search::iterator  	  Neighbour_iterator;

			typename Kernel::Compute_squared_distance_2 squared_distance;


			// Extra.
			using Log = CGAL::Level_of_detail::Mylog;

			using Cell_id  = std::pair<int, int>;
			using Grid_map = std::map<Cell_id, std::vector<int> >;

			using Grid_cell_iterator = Grid_map::const_iterator;


			// Constructor.
			Level_of_detail_grid_simplify() : 
			m_grid_cell_length(FT(1) / FT(100)),
			m_new_point_type(Grid_new_point_type::BARYCENTRE),
			m_save_result(true),
			m_silent(false) { }


			// Input functions.	
			void make_silent(const bool new_state) {
				m_silent = new_state;
			}

			void set_grid_cell_length(const FT new_cell_length) {

				assert(new_cell_length > FT(0));
				m_grid_cell_length = new_cell_length;
			}

			void set_new_point_type(const Grid_new_point_type new_point_type) {
				m_new_point_type = new_point_type;
			}

			void save_result(const bool new_state) {
				m_save_result = new_state;
			}


			// Main.
			int process(Boundary_data &, Projected_points &boundary_clutter_projected) const {

				int number_of_removed_points = 0;
				if (boundary_clutter_projected.empty()) {
					
					std::cerr << "WARNING: Grid simplify empty input!" << std::endl;
					return number_of_removed_points;
				}

				Grid_map grid_map;
				create_grid_map(grid_map, boundary_clutter_projected);

				Projected_points cleaned_points;
				simplify_clutter_using_grid(cleaned_points, grid_map, boundary_clutter_projected);

				number_of_removed_points = static_cast<int>(boundary_clutter_projected.size() - cleaned_points.size());
				boundary_clutter_projected = cleaned_points;

				assert(number_of_removed_points >= 0);
				if (m_save_result && !m_silent) {
					
					Log log; 
					log.export_projected_points_as_xyz("tmp" + std::string(PSR) + "lod_0_1" + std::string(PSR) + "grid_simplify_result", boundary_clutter_projected, "unused path");	
				}
				
				return number_of_removed_points;
			}


		// Fields.
		private:
			FT 			   		m_grid_cell_length;
			Grid_new_point_type m_new_point_type;
			bool 			    m_save_result;
			bool 				m_silent;


		// Grid simplify.
		private:
			void create_grid_map(Grid_map &grid_map, const Projected_points &boundary_clutter_projected) const {

				assert(grid_map.empty());
				assert(!boundary_clutter_projected.empty());

				Cell_id cell_id;
				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit) {
					
					const Projected_point &projected = *pit;
					get_cell_id(cell_id, projected);

					grid_map[cell_id].push_back(projected.first);
				}
			}

			void get_cell_id(Cell_id &cell_id, const Projected_point &projected) const {

				const Point_2 &point = projected.second;

				const int id_x = get_id_value(point.x());
				const int id_y = get_id_value(point.y());

				cell_id = std::make_pair(id_x, id_y);
			}

			int get_id_value(const FT value) const {

				assert(m_grid_cell_length > FT(0));
				const int id = static_cast<int>(CGAL::to_double(value / m_grid_cell_length));

				if (value >= 0) return id;
				return id - 1;
			}

			void simplify_clutter_using_grid(Projected_points &cleaned_points, const Grid_map &grid_map, const Projected_points &boundary_clutter_projected) const {

				// std::cout << "" + std::string(PN) + "Warning: grid simplify, take care when using these simplified points to access any original input data!" + std::string(PN) + "";

				assert(cleaned_points.empty());
				assert(!grid_map.empty() && !boundary_clutter_projected.empty());

				int point_index = 0;
				for (Grid_cell_iterator git = grid_map.begin(); git != grid_map.end(); ++git, ++point_index) {
					
					const auto &grid_cell = *git;
					create_new_point(cleaned_points, point_index, grid_cell, boundary_clutter_projected);
				}		
			}

			template<class Grid_cell>
			void create_new_point(Projected_points &cleaned_points, 
				const int glob_point_index, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				switch (m_new_point_type) {

					case Grid_new_point_type::CENTROID:
						create_new_centroid(cleaned_points, glob_point_index, grid_cell);
						break;

					case Grid_new_point_type::BARYCENTRE:
						create_new_barycentre(cleaned_points, glob_point_index, grid_cell, boundary_clutter_projected);
						break;

					case Grid_new_point_type::CLOSEST:
						create_new_closest(cleaned_points, grid_cell, boundary_clutter_projected);
						break;

					default:
						assert(!"Wrong new point type!");
						break;
				}
			}

			template<class Grid_cell>
			void create_new_centroid(Projected_points &cleaned_points, const int glob_point_index, const Grid_cell &grid_cell) const {

				const Cell_id &cell_id = grid_cell.first;

				const FT half_cell_length = m_grid_cell_length / FT(2);

				const FT x = static_cast<FT>(cell_id.first)  * m_grid_cell_length + half_cell_length;
				const FT y = static_cast<FT>(cell_id.second) * m_grid_cell_length + half_cell_length;

				cleaned_points[glob_point_index] = Point_2(x, y);
			}

			template<class Grid_cell>
			void create_new_barycentre(Projected_points &cleaned_points, const int glob_point_index, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				Point_2 new_point;
				create_new_barycentre_point(new_point, grid_cell, boundary_clutter_projected);

				cleaned_points[glob_point_index] = new_point;
			}

			template<class Grid_cell>
			void create_new_barycentre_point(Point_2 &new_point, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				const std::vector<int> &point_idxs = grid_cell.second;
				const size_t num_point_idxs = point_idxs.size();

				assert(num_point_idxs != 0);

				FT x = FT(0), y = FT(0);
				for (size_t i = 0; i < num_point_idxs; ++i) {

					const Point_2 &old_point = boundary_clutter_projected.at(point_idxs[i]);

					x += old_point.x();
					y += old_point.y();
				}

				x /= static_cast<FT>(num_point_idxs);
				y /= static_cast<FT>(num_point_idxs);

				new_point = Point_2(x, y);
			}

			template<class Grid_cell>
			void create_new_closest(Projected_points &cleaned_points, const Grid_cell &grid_cell, const Projected_points &boundary_clutter_projected) const {

				Point_2 barycentre;
				create_new_barycentre_point(barycentre, grid_cell, boundary_clutter_projected);

				const std::vector<int> &point_idxs = grid_cell.second;
				const size_t num_point_idxs = point_idxs.size();

				assert(num_point_idxs != 0);
				int new_point_index = -1;

				FT min_dist = FT(1000000000);
				for (size_t i = 0; i < num_point_idxs; ++i) {
					
					const Point_2 &old_point = boundary_clutter_projected.at(point_idxs[i]);
					const FT dist = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance(old_point, barycentre))));

					if (dist < min_dist) {
						
						min_dist = dist;
						new_point_index = i;
					}
				}

				assert(new_point_index >= 0 && new_point_index < static_cast<int>(num_point_idxs));
				cleaned_points[point_idxs[new_point_index]] = boundary_clutter_projected.at(point_idxs[new_point_index]);
			}

			// WARNING: DOES NOT WORK IF LATER I WILL NEED TO ACCESS ANY DATA FROM THE ORIGINAL INPUT THAT I HAD BEFORE
			// SIMPLIFICATION LIKE NORMALS OF INPUT POINTS E.G. SINCE THIS LOCAL MAPPING IS NOT VALID FOR THE ORIGINAL INPUT ANYMORE.
			void set_new_boundary_data(Boundary_data &boundary_clutter, const Projected_points &boundary_clutter_projected) const {

				boundary_clutter.clear();

				Boundary_data new_data;
				assert(boundary_clutter.empty() && !boundary_clutter_projected.empty());

				size_t count = 0;
				std::vector<int> idxs(boundary_clutter_projected.size());

				for (Point_iterator pit = boundary_clutter_projected.begin(); pit != boundary_clutter_projected.end(); ++pit, ++count) {	
					const Projected_point &projected = *pit;
					
					assert(count < idxs.size());
					idxs[count] = projected.first;
				}
				
				assert(count == idxs.size());

				new_data[0] = idxs;
				boundary_clutter = new_data;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_GRID_SIMPLIFY_H