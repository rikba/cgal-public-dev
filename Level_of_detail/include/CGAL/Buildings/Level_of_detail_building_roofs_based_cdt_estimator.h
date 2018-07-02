#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/intersections.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_based_cdt_estimator {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using Intersect = typename Kernel::Intersect_3;

            using FT      = typename Kernel::FT;
            using Line_3  = typename Kernel::Line_3;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Building_iterator = typename Buildings::iterator;

            using Roof           = typename Building::Roof;
            using Roofs          = typename Building::Roofs;
            using Envelope_input = typename Building::Data_triangles;
            using Planes         = typename Building::Planes;

            using Associated_planes = typename Roof::Associated_planes; 
            using Boundary          = typename Roof::Roof_boundary;

            Level_of_detail_building_roofs_based_cdt_estimator(const FT ground_height, Buildings &buildings) :
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_big_value(FT(100000000000000)) 
            { }

            void estimate() const {
                
                // implement!
            }

        private:
            const FT m_ground_height;
            Buildings &m_buildings;
            const FT m_big_value;
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H