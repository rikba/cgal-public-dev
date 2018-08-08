#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_creator {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            
            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            Level_of_detail_building_roofs_creator(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height)
            { }

            void create_roofs() {

                if (m_buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid && building.shapes.size() != 0 && building.jp_polygons.size() != 0)
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;

            void process_building(Building &building) const {
                
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H