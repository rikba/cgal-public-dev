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
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Building_iterator = typename Buildings::iterator;

            using Roof           = typename Building::Roof;
            using Roofs          = typename Building::Roofs;
            using Envelope_input = typename Building::Data_triangles;
            using Planes         = typename Building::Planes;
            using CDT            = typename Building::CDT;

            using Roof_boundary = typename Roof::Roof_boundary;

            using Face_handle       = typename CDT::Face_handle;
            using Vertex_handle     = typename CDT::Vertex_handle;
            using Vertices_iterator = typename CDT::Finite_vertices_iterator;
            using Face_circulator   = typename CDT::Face_circulator;

            using Associated_planes = typename Roof::Associated_planes; 
            using Boundary          = typename Roof::Roof_boundary;

            Level_of_detail_building_roofs_based_cdt_estimator(const FT ground_height, Buildings &buildings) :
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_big_value(FT(100000000000000)) 
            { }

            void estimate() {
                for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    
                    Building &building = bit->second;
                    handle_building(building);
                }
            }

        private:
            const FT m_ground_height;
            Buildings &m_buildings;
            const FT m_big_value;

            void handle_building(Building &building) const {

                const CDT &cdt = building.cdt;
                for (Vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
                    
                    const Vertex_handle &vh = static_cast<Vertex_handle>(vit);
                    handle_vertex(vh, building);
                }
            }

            void handle_vertex(const Vertex_handle &vh, Building &building) const {

                const Point_2 &query = vh->point();
                const Line_3 line = Line_3(Point_3(query.x(), query.y(), FT(0)), Point_3(query.x(), query.y(), FT(10)));

                const CDT &cdt = building.cdt;
                Roofs &roofs   = building.roofs;

                Face_circulator face_circulator = cdt.incident_faces(vh);
                Face_circulator end = face_circulator;

                if (face_circulator.is_empty()) return;

                std::vector<FT> zs;
                do {
                    Face_handle face = static_cast<Face_handle>(face_circulator);
                    const int roof_index = face->info().roof_index;

                    if (cdt.is_infinite(face) || roof_index < 0) {
                        
                        ++face_circulator; 
                        continue;
                    }
                    
                    const Roof &roof = building.roofs[roof_index];
                    const Associated_planes &planes = roof.associated_planes;

                    FT z;
                    const size_t index = planes[0];

                    if (!roof.is_plane_index) {
                        if (building.envelope_input[index].second.is_vertical) {

                            ++face_circulator;
                            continue;
                        }

                        const Point_3 &p1 = building.envelope_input[index].first.vertex(0);
                        const Point_3 &p2 = building.envelope_input[index].first.vertex(1);
                        const Point_3 &p3 = building.envelope_input[index].first.vertex(2);

                        const Plane_3 plane = Plane_3(p1, p2, p3);
                        z = intersect_line_with_plane(line, plane);

                    } else {

                        const Plane_3 &plane = building.planes[index];
                        z = intersect_line_with_plane(line, plane);
                    }

                    zs.push_back(z);
                    ++face_circulator;

                } while (face_circulator != end);
                
                // if (zs.size() == 0) return;

                FT final_z = FT(0);
                for (size_t i = 0; i < zs.size(); ++i) final_z += zs[i];
                final_z /= static_cast<FT>(zs.size());

                face_circulator = end;
                do {
                    Face_handle face = static_cast<Face_handle>(face_circulator);
                    const int roof_index = face->info().roof_index;

                    if (cdt.is_infinite(face) || roof_index < 0) {
                        
                        ++face_circulator; 
                        continue;
                    }
                    
                    Roof &roof = building.roofs[roof_index];
                    Roof_boundary &boundary = roof.boundary;

                    const Associated_planes &planes = roof.associated_planes;
                    const size_t index = planes[0];

                    if (building.envelope_input[index].second.is_vertical) {

                        ++face_circulator;
                        continue;
                    }

                    Point_3 &p1 = boundary[0];
                    Point_3 &p2 = boundary[1];
                    Point_3 &p3 = boundary[2];

                    const FT eps = FT(1) / FT(100);

                    if (CGAL::abs(p1.x() - query.x()) < eps && CGAL::abs(p1.y() - query.y()) < eps) {
                        p1 = Point_3(p1.x(), p1.y(), final_z);

                        ++face_circulator; 
                        continue;
                    }

                    if (CGAL::abs(p2.x() - query.x()) < eps && CGAL::abs(p2.y() - query.y()) < eps) {
                        p2 = Point_3(p2.x(), p2.y(), final_z);
                        
                        ++face_circulator; 
                        continue;
                    }

                    if (CGAL::abs(p3.x() - query.x()) < eps && CGAL::abs(p3.y() - query.y()) < eps) {
                        p3 = Point_3(p3.x(), p3.y(), final_z);
                        
                        ++face_circulator; 
                        continue;
                    }

                    ++face_circulator;
                } while (face_circulator != end);
            }

            FT intersect_line_with_plane(const Line_3 &line, const Plane_3 &plane) const {

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                const Point_3 r = boost::get<Point_3>(*result);

				return r.z();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H