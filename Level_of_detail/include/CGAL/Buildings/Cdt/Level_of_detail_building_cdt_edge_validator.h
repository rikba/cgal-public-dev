#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_EDGE_VALIDATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_EDGE_VALIDATOR_H

// CGAL includes.
#include <CGAL/utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding>
		class Level_of_detail_building_cdt_edge_validator {
            
        public:
            typedef InputKernel   Kernel;
            typedef InputBuilding Building;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;

            using CDT         = typename Building::CDT;
            using Face_handle = typename CDT::Face_handle;

            Level_of_detail_building_cdt_edge_validator(Building &building) :
            m_building(building)
            { }

            void validate() const {
                
            }

        private:
            Building &m_building;
        };
        
    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_EDGE_VALIDATOR_H