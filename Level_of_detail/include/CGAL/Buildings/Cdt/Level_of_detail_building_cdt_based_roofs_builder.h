#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_BASED_ROOFS_BUILDER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_BASED_ROOFS_BUILDER_H

// New CGAL includes.
#include <CGAL/Buildings/Associaters/Level_of_detail_building_partition_vote_based_plane_associater.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_cdt_based_roofs_builder {

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

            using Face_handle    = typename CDT::Face_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;

            using Roof  = typename Building::Roof;
            using Roofs = typename Building::Roofs;

            using Plane_associater = CGAL::LOD::Level_of_detail_building_partition_vote_based_plane_associater<Kernel, Input, Building>;

            Level_of_detail_building_cdt_based_roofs_builder(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_ground_height(ground_height),
            m_buildings(buildings)
            { }

            void build() {

                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    if (!building.is_valid) continue;
                    build_roofs(building);
                }
            }

        private:
            const Input &m_input;
            const FT m_ground_height;
            Buildings &m_buildings;

            void build_roofs(Building &building) const {

                add_roof_faces(building);
                add_associated_planes(building);
            }

            void add_roof_faces(Building &building) const {
                
                const CDT &cdt = building.cdt;
                const FT  z    = building.height + m_ground_height;
                
                Roofs &roofs = building.roofs;
                roofs.clear();

                int roof_index = 0; Roof roof;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

                    roof.boundary.clear();
                    Face_handle fh = static_cast<Face_handle>(fit);

                    if (!fh->info().is_valid) continue;

                    const Point_2 &p1 = fh->vertex(0)->point();
                    const Point_2 &p2 = fh->vertex(1)->point();
                    const Point_2 &p3 = fh->vertex(2)->point();

                    roof.boundary.push_back(Point_3(p1.x(), p1.y(), z));
                    roof.boundary.push_back(Point_3(p2.x(), p2.y(), z));
                    roof.boundary.push_back(Point_3(p3.x(), p3.y(), z));

                    roof.index = roof_index; ++roof_index;
                    fh->info().roof_index = roof.index;

                    roofs.push_back(roof);
                }
            }

            void add_associated_planes(Building &building) const {
                
                const FT reference_height = building.roofs_min_height + FT(1) / FT(2);
                Plane_associater plane_associater(m_input, building, reference_height);

                for (size_t i = 0; i < building.roofs.size(); ++i) {
                    
                    building.roofs[i].associated_planes.clear();
                    plane_associater.find_associated_planes(i, building.roofs[i].is_plane_index, building.roofs[i].associated_planes);
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_BASED_ROOFS_BUILDER_H