#ifndef CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding>
		class Level_of_detail_polyhedron_in_out_estimator {

		public:
            using Kernel   = InputKernel;
			using Input    = InputData;
			using Building = InputBuilding;

			using FT 	     = typename Kernel::FT;
			using Line_2     = typename Kernel::Line_2;
			using Point_2    = typename Kernel::Point_2;
			using Point_3    = typename Kernel::Point_3;
			using Triangle_2 = typename Kernel::Triangle_2;

            using Index       = typename Building::Index;
            using Indices     = typename Building::Indices;
			using Polyhedron  = typename Building::Polyhedron;
			using Polyhedrons = typename Building::Polyhedrons;
			using Floor_faces = typename Building::Floor_faces;

			using Vertices = typename Polyhedron::Vertices;

			using Pair = CGAL::cpp11::array<FT, 2>;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

			Level_of_detail_polyhedron_in_out_estimator(const Input &input, const FT ground_height) : 
			m_input(input),
			m_ground_height(ground_height),
			m_big_value(FT(100000000000000)),
			m_distance_tolerance(FT(2)),
			m_bc_tolerance_top(FT(6) / FT(5)),
            m_bc_tolerance_bottom(-FT(1) / FT(5)) { }

			void estimate(Building &building) const {
				compute_building_maximum_height(building);

				Polyhedrons &polyhedrons = building.polyhedrons;
                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    
					Polyhedron &polyhedron = polyhedrons[i];
					polyhedron.out = estimate_out_value(building, polyhedron);

					CGAL_precondition(polyhedron.out >= FT(0) && polyhedron.out <= FT(1));
					polyhedron.in = FT(1) - polyhedron.out;
                }
			}

		private:
			const Input &m_input;

			const FT m_ground_height;
			const FT m_big_value;
			const FT m_distance_tolerance;
			const FT m_bc_tolerance_top;
            const FT m_bc_tolerance_bottom;

            void compute_building_maximum_height(Building &building) const {

                const Indices &interior_indices = building.interior_indices;
                CGAL_precondition(interior_indices.size() > 0);

                FT max_height = -m_big_value;
				for (size_t i = 0; i < interior_indices.size(); ++i) {
                    
                    const Index point_index = interior_indices[i];
                    const Point_3 &p = m_input.point(point_index);

                    max_height = CGAL::max(max_height, p.z());
                }

                building.max_height = max_height;
                CGAL_postcondition(max_height >= FT(0));
            }

			FT estimate_out_value(const Building &building, const Polyhedron &polyhedron) const {

				FT in  = FT(0);
				FT out = FT(0);

				Point_3 barycentre;
                compute_polyhedron_barycentre(polyhedron, barycentre);

                if (is_above_building_max_height(barycentre, building.max_height)) return FT(1);
                if (is_below_ground(barycentre, m_ground_height)) 				   return FT(1);
                if (is_out_of_building(barycentre, building)) 					   return FT(1);

                // if (has_vertices_outside(polyhedron, building)) out += FT(1);
                // if (is_statistically_invalid(polyhedron, building)) out += FT(1);
                
                return FT(0);
			}

			void compute_polyhedron_barycentre(const Polyhedron &polyhedron, Point_3 &barycentre) const {

                const Vertices &vertices  = polyhedron.vertices;
                const size_t num_vertices = vertices.size();

                CGAL_precondition(num_vertices > 0);

                FT x = FT(0), y = FT(0), z = FT(0);
                for (size_t i = 0; i < num_vertices; ++i) {

                    x += vertices[i].x();
                    y += vertices[i].y();
                    z += vertices[i].z();
                }

                x /= static_cast<FT>(num_vertices);
                y /= static_cast<FT>(num_vertices);
                z /= static_cast<FT>(num_vertices);

                barycentre = Point_3(x, y, z);
            }

			bool is_above_building_max_height(const Point_3 &query, const FT building_max_height) const {
                return query.z() > building_max_height;
            }

			bool is_below_ground(const Point_3 &query, const FT ground_height) const {
                return query.z() < ground_height;
            }

			bool is_out_of_building(const Point_3 &query, const Building &building) const {
                
                const Point_2 p = Point_2(query.x(), query.y());

                const Floor_faces &floor_faces = building.faces;
                for (size_t i = 0; i < floor_faces.size(); ++i) {
                        
                    const Point_2 &p1 = floor_faces[i]->vertex(0)->point();
                    const Point_2 &p2 = floor_faces[i]->vertex(1)->point();
                    const Point_2 &p3 = floor_faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (is_within_triangle(p, triangle)) return false;
                }
                return true;
            }

			bool is_within_triangle(const Point_2 &query, const Triangle_2 &triangle) const {
                
                if (triangle.has_on_bounded_side(query) || triangle.has_on_boundary(query)) 
                    return true;
                
                for (size_t i = 0; i < 3; ++i) {
                    const size_t ip = (i + 1) % 3;

                    const Point_2 &p1 = triangle.vertex(i);
                    const Point_2 &p2 = triangle.vertex(ip);

                    const Line_2 line = Line_2(p1, p2);

                    const Point_2 projected   = line.projection(query);
                    const FT squared_distance = squared_distance_2(query, projected);

                    const Pair pair = BC::compute_segment_coordinates_2(p1, p2, projected, Kernel());

                    const FT squared_tolerance = m_distance_tolerance * m_distance_tolerance;
                    
                    const FT epst = m_bc_tolerance_top;
                    const FT epsb = m_bc_tolerance_bottom;

                    if (pair[0] > epsb && pair[1] > epsb && pair[0] < epst && pair[1] < epst && squared_distance < squared_tolerance) return true;
                }
                return false;
            }
		};

	} // CGAL

} // LOD

#endif // CGAL_LEVEL_OF_DETAIL_POLYHEDRON_IN_OUT_ESTIMATOR_H