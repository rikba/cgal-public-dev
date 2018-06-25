#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>

// New CGAL includes.
#include <CGAL/Buildings/Level_of_detail_building_roof_face_validator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_based_cdt_cleaner {

        public:
            using Kernel    = InputKernel;
            using Input     = InputData;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;

            using CDT                = typename Building::CDT;
            using Buildings_iterator = typename Buildings::iterator;

            using Edge           = typename CDT::Edge;
            using Face_handle    = typename CDT::Face_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;

            using Label        = int;
            using Labels       = std::vector<Label>;
            using Label_pair   = std::pair<Face_handle, Labels>;
            using Sorted_faces = std::map<Face_handle, Labels>;
            using Comparator   = std::function<bool(Label_pair, Label_pair)>;
            using Label_set    = std::set<Label_pair, Comparator>;

            using Index   = typename Building::Index;
			using Indices = typename Building::Indices;
            using Shapes  = typename Building::Shapes;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;
            
            using Locate_type         = typename CDT::Locate_type;
            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;

            Level_of_detail_building_roofs_based_cdt_cleaner(const Input &input, const FT ground_height, Buildings &buildings) : 
            m_input(input),
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_thin_face_max_size(FT(1) / FT(2)),
            m_max_percentage(FT(95)) 
            { }

            void clean() {
                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    Building &building = (*bit).second;
                    
                    clean_cdt(building);
                    m_roof_face_validator.mark_wrong_faces(building);
                }
            }

        private:
            const Input &m_input;
            const FT m_ground_height;

            Buildings &m_buildings;

            const FT m_thin_face_max_size;
            const FT m_max_percentage;

            Roof_face_validator m_roof_face_validator;

            void clean_cdt(Building &building) const {

                CDT &cdt = building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    Face_handle fh = static_cast<Face_handle>(fit);

                    if (!fh->info().is_valid) continue;
                    handle_face(fh, cdt);
                }
            }

            void handle_face(Face_handle &fh, CDT &cdt) const {


            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H