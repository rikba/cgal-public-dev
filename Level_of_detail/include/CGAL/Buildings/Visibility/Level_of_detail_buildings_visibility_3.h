#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_buildings_visibility_3 {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_2 = typename Kernel::Triangle_2;
            
            using Buildings_iterator = typename Buildings::iterator;

            using Index       = typename Building::Index;
            using Indices     = typename Building::Indices;
            using Floor_faces = typename Building::Floor_faces;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Vertices = typename Polyhedron::Vertices;

            Level_of_detail_buildings_visibility_3(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_big_value(FT(100000000000000))
            { }

            void apply() {

                if (m_buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) {

                        compute_building_maximum_height(building);
                        process_building(building);
                    }
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;
            const FT m_big_value;

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

            void process_building(Building &building) const {

                Polyhedrons &polyhedrons = building.polyhedrons;
                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    
                    Polyhedron &polyhedron = polyhedrons[i];
                    polyhedron.is_valid = is_valid_polyhedron(building, polyhedron);
                }
            }

            bool is_valid_polyhedron(const Building &building, const Polyhedron &polyhedron) const {

                Point_3 b;
                compute_polyhedron_barycentre(polyhedron, b);

                if (is_above_building(b, building.max_height)) return false;
                if (is_below_ground(b, m_ground_height)) return false;
                if (is_outside_building(b, building)) return false;

                return true;
            }

            void compute_polyhedron_barycentre(const Polyhedron &polyhedron, Point_3 &b) const {

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

                b = Point_3(x, y, z);
            }

            bool is_above_building(const Point_3 &query, const FT building_max_height) const {
                return query.z() > building_max_height;
            }

            bool is_below_ground(const Point_3 &query, const FT ground_height) const {
                return query.z() < ground_height;
            }

            bool is_outside_building(const Point_3 &query, const Building &building) const {
                
                const Point_2 p = Point_2(query.x(), query.y());

                const Floor_faces &floor_faces = building.faces;
                for (size_t i = 0; i < floor_faces.size(); ++i) {
                        
                    const Point_2 &p1 = floor_faces[i]->vertex(0)->point();
                    const Point_2 &p2 = floor_faces[i]->vertex(1)->point();
                    const Point_2 &p3 = floor_faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (triangle.has_on_bounded_side(p) || triangle.has_on_boundary(p)) return false;
                }
                return true;
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H