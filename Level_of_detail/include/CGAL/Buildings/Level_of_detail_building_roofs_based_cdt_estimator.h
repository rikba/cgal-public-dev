#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H

// STL includes.
#include <map>
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

            using FT         = typename Kernel::FT;
            using Line_3     = typename Kernel::Line_3;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Triangle_3 = typename Kernel::Triangle_3;

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

            using Roof_indices  = std::vector<int>;
            using Heights       = std::map<int, FT>;
            using Final_height  = std::pair<FT, Roof_indices>;
            using Final_heights = std::vector<Final_height>;

            using Planes_iterator        = typename Planes::const_iterator;
            using Heights_iterator       = typename Heights::const_iterator; 
            using Final_heights_iterator = typename Final_heights::const_iterator;

            Level_of_detail_building_roofs_based_cdt_estimator(const FT ground_height, Buildings &buildings) :
            m_tolerance(FT(1) / FT(1000000)),
            m_ground_height(ground_height),
            m_buildings(buildings)
            { }

            void estimate() {
                for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    
                    Building &building = bit->second;
                    handle_building(building);
                }
            }

        private:
            const FT m_tolerance;
            const FT m_ground_height;

            Buildings &m_buildings;
            
            void handle_building(Building &building) const {

                const CDT &cdt = building.cdt;
                for (Vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
                    
                    const Vertex_handle &vh = static_cast<Vertex_handle>(vit);
                    handle_vertex(vh, building);
                }
            }

            void handle_vertex(const Vertex_handle &vh, Building &building) const {

                Line_3 line;
                get_line(vh, line);

                Planes planes; Roof_indices roof_indices;
                get_planes_with_indices(vh, building, planes, roof_indices);

                CGAL_precondition(roof_indices.size() == planes.size());
                if (planes.size() == 0) return;

                Heights heights;
                compute_heights(line, planes, building.height, roof_indices, heights);

                if (heights.size() == 0) return;

                Final_heights final_heights;
                compute_final_heights(heights, final_heights);

                if (final_heights.size() == 0) return;

                update_roof_heights(vh, final_heights, building);
            }

            void get_line(const Vertex_handle &vh, Line_3 &line) const {
                
                const Point_2 &point = vh->point();
                line = Line_3(Point_3(point.x(), point.y(), FT(0)), Point_3(point.x(), point.y(), FT(10)));
            }

            void get_planes_with_indices(const Vertex_handle &vh, const Building &building, Planes &planes, Roof_indices &roof_indices) const {

                planes.clear();
                roof_indices.clear();

                const CDT &cdt     = building.cdt;
                const Roofs &roofs = building.roofs;

                Face_circulator face_circulator = cdt.incident_faces(vh);
                Face_circulator end = face_circulator;

                if (face_circulator.is_empty()) return;

                do {
                    const Face_handle &face = static_cast<Face_handle>(face_circulator);
                    const int roof_index = face->info().roof_index;

                    if (cdt.is_infinite(face) || roof_index < 0) {
                        ++face_circulator; 

                        if (face_circulator == end) break;
                        else continue;
                    }

                    const Roof &roof = roofs[roof_index];

                    if (!roof.is_valid) {
                        ++face_circulator; 
                        
                        if (face_circulator == end) break;
                        else continue;
                    }

                    const Associated_planes &plane_indices = roof.associated_planes;

                    if (plane_indices.size() < 1) {
                        ++face_circulator;

                        if (face_circulator == end) break;
                        else continue;
                    }

                    const size_t index = plane_indices[0];

                    if (!roof.is_plane_index) {
                        const Envelope_input &triangles = building.envelope_input;

                        if (triangles[index].second.is_vertical) {
                            ++face_circulator;
                            
                            if (face_circulator == end) break;
                            else continue;
                        }

                        const Triangle_3 &triangle = triangles[index].first;

                        const Point_3 &p1 = triangle.vertex(0);
                        const Point_3 &p2 = triangle.vertex(1);
                        const Point_3 &p3 = triangle.vertex(2);

                        const Plane_3 plane = Plane_3(p1, p2, p3);
                        planes.push_back(plane);

                    } else {

                        const Planes &associated_planes = building.planes; 
                        const Plane_3 &plane = associated_planes[index];
                        planes.push_back(plane);
                    }

                    roof_indices.push_back(roof_index);
                    ++face_circulator;

                } while (face_circulator != end);
            }

            void compute_heights(const Line_3 &line, const Planes &planes, const FT building_height, const Roof_indices &roof_indices, Heights &heights) const {
                
                CGAL_precondition(roof_indices.size() == planes.size());
                heights.clear();

                size_t i = 0;
                for (Planes_iterator pit = planes.begin(); pit != planes.end(); ++pit, ++i) {
                    const Plane_3 &plane = *pit;

                    const FT height = intersect_line_with_plane(line, plane);
                    if (height > m_ground_height + building_height / FT(2) && height < m_ground_height + building_height * FT(2)) heights[roof_indices[i]] = height;
                }
            }

            FT intersect_line_with_plane(const Line_3 &line, const Plane_3 &plane) const {

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                const Point_3 point = boost::get<Point_3>(*result);

				return point.z();
            }

            void compute_final_heights(const Heights &heights, Final_heights &final_heights) const {
                
                final_heights.clear();
                const FT final_height = compute_average_height(heights);
                
                Roof_indices final_indices;
                add_all_indices(heights, final_indices);

                final_heights.push_back(std::make_pair(final_height, final_indices));
            }

            FT compute_average_height(const Heights &heights) const {
                
                FT avg_height = FT(0);
                for (Heights_iterator hit = heights.begin(); hit != heights.end(); ++hit) {
                    
                    const FT height = hit->second;
                    avg_height += height;
                }
                
                avg_height /= static_cast<FT>(heights.size());
                return avg_height;
            }

            void add_all_indices(const Heights &heights, Roof_indices &roof_indices) const {
                
                roof_indices.clear();
                for (Heights_iterator hit = heights.begin(); hit != heights.end(); ++hit) {
                    
                    const int roof_index = hit->first;
                    roof_indices.push_back(roof_index);
                }
            }

            void update_roof_heights(const Vertex_handle &vh, const Final_heights &final_heights, Building &building) const {

                Roofs &roofs = building.roofs;
                for (Final_heights_iterator fh_it = final_heights.begin(); fh_it != final_heights.end(); ++fh_it) {
                    
                    const FT final_height = fh_it->first; 
                    const Roof_indices &roof_indices = fh_it->second;

                    for (size_t i = 0; i < roof_indices.size(); ++i)
                        update_roof_height(vh->point(), final_height, roofs[roof_indices[i]]);
                }
            }

            void update_roof_height(const Point_2 &ref, const FT final_height, Roof &roof) const {
                Roof_boundary &roof_boundary = roof.boundary;
                
                for (size_t i = 0; i < roof_boundary.size(); ++i) {
                    Point_3 &point = roof_boundary[i];

                    if (are_equal(point, ref)) {

                        update_point_height(final_height, point);
                        return;
                    }
                }
            }

            inline bool are_equal(const Point_3 &p, const Point_2 &ref) const {
                return CGAL::abs(p.x() - ref.x()) < m_tolerance && CGAL::abs(p.y() - ref.y()) < m_tolerance;
            }

            inline void update_point_height(const FT final_height, Point_3 &point) const {
                point = Point_3(point.x(), point.y(), final_height);
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_ESTIMATOR_H