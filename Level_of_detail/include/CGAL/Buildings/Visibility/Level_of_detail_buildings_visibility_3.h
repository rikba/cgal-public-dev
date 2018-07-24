#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Buildings/Roofs/Roof/Level_of_detail_building_roof_face_validator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_buildings_visibility_3 {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT = typename Kernel::FT;
            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Vertex   = typename Polyhedron::Vertex;
            using Vertices = typename Polyhedron::Vertices;

            using Facet  = typename Polyhedron::Facet;
            using Facets = typename Polyhedron::Facets;

            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;
            using Boundary = typename Roof_face_validator::Boundary;

            Level_of_detail_buildings_visibility_3(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_ground_tolerance(FT(1) / FT(10))
            { }

            void apply() {

                if (m_buildings.size() == 0)
                    return; 
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;
            const FT m_ground_tolerance;

            Roof_face_validator m_roof_face_validator;

            void process_building(Building &building) const {

                Polyhedrons &polyhedrons = building.polyhedrons;
                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    
                    Polyhedron &polyhedron = polyhedrons[i];
                    check_polyhedron_validity(building, polyhedron);
                }
            }

            void check_polyhedron_validity(const Building &building, Polyhedron &polyhedron) const {
                
                if (!is_valid_polyhedron(building, polyhedron)) {
                    
                    polyhedron.is_valid = false;
                    return;
                }

                polyhedron.is_valid = true;
            }

            bool is_valid_polyhedron(const Building &building, Polyhedron &polyhedron) const {

                const Vertices &vertices = polyhedron.vertices;
                const Facets     &facets = polyhedron.facets;

                // Search for wrong exterior polyhedrons.
                for (size_t i = 0; i < facets.size(); ++i) {        
                    const Facet &facet = facets[i];

                    Boundary boundary(facet.size());
                    for (size_t j = 0; j < facet.size(); ++j) boundary[j] = vertices[facet[j]];

                    if (!m_roof_face_validator.is_valid_polyhedron_facet(building, boundary, true)) return false;
                }

                // Search for wrong interior polyhedrons.
                size_t count = 0;
                for (size_t i = 0; i < vertices.size(); ++i) {
                    const Vertex &vertex = vertices[i];

                    const FT vertex_height = vertex.z();
                    if (CGAL::abs(vertex_height - m_ground_height) > m_ground_tolerance) ++count;
                }

                return count != vertices.size();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_VISIBILITY_3_H