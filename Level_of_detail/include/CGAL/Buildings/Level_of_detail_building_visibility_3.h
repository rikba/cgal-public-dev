#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_VISIBILITY_3_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_VISIBILITY_3_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Buildings/Level_of_detail_building_roof_face_validator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputCDT, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_visibility_3 {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using CDT       = InputCDT;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT = typename Kernel::FT;
            using Building_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Vertex   = typename Polyhedron::Vertex;
            using Vertices = typename Polyhedron::Vertices;

            using Facet  = typename Polyhedron::Facet;
            using Facets = typename Polyhedron::Facets;

            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;
            using Boundary = typename Roof_face_validator::Boundary;

            Level_of_detail_building_visibility_3(const Input &input, const CDT &cdt, const FT ground_height) :
            m_input(input),
            m_cdt(cdt),
            m_ground_height(ground_height)
            { }

            void filter(Buildings &buildings) const {

                if (buildings.size() == 0) return; int count = 0;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit, ++count) {

                    Building &building = bit->second; building.index = count;
					if (building.is_valid) process_building(building);
                }
            }

        private:
            const Input &m_input;
            const CDT   &m_cdt;

            const FT m_ground_height;

            Roof_face_validator m_roof_face_validator;

            void process_building(Building &building) const {

                Polyhedrons &polyhedrons = building.polyhedrons;
                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    
                    Polyhedron &polyhedron = polyhedrons[i];
                    check_polyhedron_validity(building, polyhedron);
                }
            }

            void check_polyhedron_validity(const Building &building, Polyhedron &polyhedron) const {
                
                if (!is_vertex_criterion(building, polyhedron)) {
                    
                    polyhedron.is_valid = false;
                    return;
                }

                polyhedron.is_valid = true;
            }

            bool is_vertex_criterion(const Building &building, Polyhedron &polyhedron) const {
                
                const Vertices &vertices = polyhedron.vertices;
                const Facets     &facets = polyhedron.facets;

                for (size_t i = 0; i < facets.size(); ++i) {        
                    const Facet &facet = facets[i];

                    Boundary boundary(facet.size());
                    for (size_t j = 0; j < facet.size(); ++j) boundary[j] = vertices[facet[j]];

                    if (!m_roof_face_validator.is_valid_polyhedron_facet(building, boundary, true)) return false;
                }

                return true;
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_VISIBILITY_3_H