#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CREATOR_H

// STL includes.
#include <vector>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_based_cdt_creator {

        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using CDT                = typename Building::CDT;
            using Buildings_iterator = typename Buildings::iterator;

            using Vertex_handle = typename CDT::Vertex_handle;
            using Roof          = typename Building::Roof;
            using Roofs         = typename Building::Roofs;
            using Roof_boundary = typename Roof::Roof_boundary;

            using Roofs_iterator         = typename Roofs::const_iterator;
            using Roof_boundary_iterator = typename Roof_boundary::const_iterator;

            using Vertex_handles = std::vector< std::vector<Vertex_handle> >;

            Level_of_detail_building_roofs_based_cdt_creator(Buildings &buildings) : 
            m_buildings(buildings) 
            { }

            void create() {
                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    
                    Building &building = (*bit).second;
                    create_cdt(building);
                }
            }

        private:
            Buildings &m_buildings;

            void create_cdt(Building &building) const {

                CDT &cdt = building.cdt;
                cdt.clear();

                const Roofs &roofs = building.roofs;
                Vertex_handles vhs(roofs.size());

                // Insert points.
                size_t i = 0;
                for (Roofs_iterator rit = roofs.begin(); rit != roofs.end(); ++rit, ++i) {
                    const Roof_boundary &roof_boundary = rit->boundary;

                    vhs[i].resize(roof_boundary.size()); size_t j = 0;
                    for (Roof_boundary_iterator rb_it = roof_boundary.begin(); rb_it != roof_boundary.end(); ++rb_it, ++j) {
                        
                        const Point_3 &point = *rb_it;
                        vhs[i][j] = cdt.insert(Point_2(point.x(), point.y()));
                    }
                }

                // Insert constraints.
                const size_t size_i = vhs.size();
                for (size_t i = 0; i < size_i; ++i){

                    const size_t size_j = vhs[i].size();
					for (size_t j = 0; j < size_j; ++j) {

						const size_t jp = (j + 1) % size_j;
						if (vhs[i][j] != vhs[i][jp]) cdt.insert_constraint(vhs[i][j], vhs[i][jp]);
					}
				}
            }
        };
        
    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CREATOR_H