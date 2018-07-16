#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_CLEANER_H

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// New CGAL includes.
#include <CGAL/Buildings/Cdt/Level_of_detail_building_cdt_based_roofs_builder.h>
#include <CGAL/Buildings/Cdt/Level_of_detail_building_cdt_edge_validator.h>
#include <CGAL/Buildings/Cdt/Level_of_detail_building_cdt_face_validator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_cdt_cleaner {

        public:
            using Kernel    = InputKernel;
            using Input     = InputData;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using CDT                = typename Building::CDT;
            using Buildings_iterator = typename Buildings::iterator;

            using Edge            = typename CDT::Edge;
            using Face_handle     = typename CDT::Face_handle;
            using Vertex_handle   = typename CDT::Vertex_handle;
            using Faces_iterator  = typename CDT::Finite_faces_iterator;
            using Edge_circulator = typename CDT::Edge_circulator;
            using Face_circulator = typename CDT::Face_circulator;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using Builder = CGAL::LOD::Level_of_detail_building_cdt_based_roofs_builder<Kernel, Input, Building, Buildings>;
            
            using Edge_validator = CGAL::LOD::Level_of_detail_building_cdt_edge_validator<Kernel, Building>;
            using Face_validator = CGAL::LOD::Level_of_detail_building_cdt_face_validator<Kernel, Building>;

            Level_of_detail_building_cdt_cleaner(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_ground_height(ground_height),
            m_buildings(buildings)
            { }

            void clean() {

                clean_cdt();
                build_data_structure();
            }

        private:
            const Input &m_input;
            const FT m_ground_height;
            Buildings &m_buildings;

            std::shared_ptr<Face_validator> m_face_validator;

            void clean_cdt() {
                
                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    if (!building.is_valid) continue;
                    clean_cdt(building);

                    m_face_validator = std::make_shared<Face_validator>(building);
                    m_face_validator->validate();
                }
            }

            void clean_cdt(Building &building) const {

            }

            void build_data_structure() {

                Builder builder(m_input, m_ground_height, m_buildings);
                builder.build();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_CLEANER_H