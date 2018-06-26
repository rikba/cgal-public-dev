#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
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

            using Edge            = typename CDT::Edge;
            using Face_handle     = typename CDT::Face_handle;
            using Vertex_handle   = typename CDT::Vertex_handle;
            using Faces_iterator  = typename CDT::Finite_faces_iterator;
            using Edge_circulator = typename CDT::Edge_circulator;

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

            using Color = CGAL::Color;
            using Log   = CGAL::LOD::Mylog;
            
            using Vertex_handles = std::vector<Vertex_handle>;

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

                bool collapse_detected = false;
                do {

                    CDT &cdt = building.cdt; int count = 0;
                    for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
                        fit->info().index = count++;

                    for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                        Face_handle fh = static_cast<Face_handle>(fit);

                        if (fh->info().is_checked && !fh->info().is_valid) continue;
                        collapse_detected = handle_face(fh, building);

                        if (collapse_detected)
                            break;
                    }

                } while (collapse_detected);
            }

            bool handle_face(const Face_handle &fh, Building &building) const {
                
                if (should_be_collapsed(fh)) return collapse_face(fh, building);
                return false;
            }

            bool should_be_collapsed(const Face_handle &fh) const {
                
                if (is_thin(fh)) return true;
                if (is_insignificant(fh)) return true;

                return false;
            }

            inline bool is_thin(const Face_handle &fh) const {
                return m_roof_face_validator.is_ghost_face(fh, m_thin_face_max_size);
            }

            bool is_insignificant(const Face_handle &fh) const {
                return false;
            }

            bool collapse_face(const Face_handle &fh, Building &building) const {

                const int neighbour_index = is_boundary_face(fh, building);
                if (neighbour_index >= 0) {
                    const Edge boundary_edge = std::make_pair(fh, neighbour_index);
                    
                    collapse_boundary_face(boundary_edge, building);
                    return true;
                }
                
                collapse_interior_face(fh, building);
                return false;
            }

            int is_boundary_face(const Face_handle &fh, const Building &building) const {
                
                for (int i = 0; i < 3; ++i) {
                    Face_handle fhn = fh->neighbor(i);
                    
                    if (fhn->info().is_checked && !fhn->info().is_valid) return i;
                    
                    if (building.cdt.is_infinite(fhn)) return i;
                    
                    if (!fhn->info().is_checked) {
                        fhn->info().is_valid = m_roof_face_validator.is_valid_face(building, fhn);
                        
                        if (!fhn->info().is_valid) 
                            return i;
                    }
                }
                return -1;
            }

            void collapse_boundary_face(const Edge &boundary_edge, Building &building) const {
                CGAL_precondition(boundary_edge.second >= 0 && boundary_edge.second < 3);

                Vertex_handle vh1, vh2;
                find_edge_vertices(boundary_edge, vh1, vh2);

                Point_2 new_point;
                find_new_edge_point(vh1, vh2, new_point);

                Vertex_handles vhs1;
                update_incident_constraints(vh1, vhs1, building);

                Vertex_handles vhs2;
                update_incident_constraints(vh2, vhs2, building);

                Vertex_handles vhs;
                merge_constraints(vhs1, vhs2, vh1, vh2, vhs);

                building.cdt.remove(vh1);
                building.cdt.remove(vh2);

                const Vertex_handle vh = building.cdt.insert(new_point);
                insert_constraints(vh, vhs, building);
            }

            void find_edge_vertices(const Edge &edge, Vertex_handle &vh1, Vertex_handle &vh2) const {
                const int vertex_index = edge.second;

                const int vi1 = (vertex_index + 2) % 3;
                const int vi2 = (vertex_index + 1) % 3;

                vh1 = edge.first->vertex(vi1);
                vh2 = edge.first->vertex(vi2);
            }

            inline void find_new_edge_point(const Vertex_handle &vh1, const Vertex_handle &vh2, Point_2 &new_point) const {
                find_mid_point(vh1, vh2, new_point);
            }

            void find_mid_point(const Vertex_handle &vh1, const Vertex_handle &vh2, Point_2 &result) const {

                const Point_2 &p1 = vh1->point();
                const Point_2 &p2 = vh2->point();

                const FT x = (p1.x() + p2.x()) / FT(2);
                const FT y = (p1.y() + p2.y()) / FT(2);

                result = Point_2(x, y);
            }

            void update_incident_constraints(const Vertex_handle &query, Vertex_handles &vhs, Building &building) const {
                
                vhs.clear();
                CDT &cdt = building.cdt;

                Edge_circulator edge_circulator = cdt.incident_edges(query);
                Edge_circulator end = edge_circulator;

                do {
                    const Edge &edge = *edge_circulator;
                    if (cdt.is_constrained(edge) && !cdt.is_infinite(edge)) {

                        Vertex_handle vh1, vh2;
                        find_edge_vertices(edge, vh1, vh2);

                        CGAL_precondition(vh1 == query || vh2 == query);
                        if (vh1 != query) vhs.push_back(vh1);
                        if (vh2 != query) vhs.push_back(vh2);

                        cdt.remove_constraint(edge.first, edge.second);
                    }

                    ++edge_circulator;
                } while (edge_circulator != end);
            }

            void merge_constraints(
            const Vertex_handles &vhs1, const Vertex_handles &vhs2, 
            const Vertex_handle   &vh1, const Vertex_handle   &vh2, Vertex_handles &result) const {
                result.clear();

                add_handles(vhs1, vh2, result);
                add_handles(vhs2, vh1, result);
            }

            void add_handles(const Vertex_handles &vhs, const Vertex_handle &bad, Vertex_handles &result) const {
                
                for (size_t i = 0; i < vhs.size(); ++i) {
                    const Vertex_handle &query = vhs[i];

                    if (query != bad && handle_can_be_included(result, query)) 
                        result.push_back(query);
                }
            }

            inline bool handle_can_be_included(const Vertex_handles &vhs, const Vertex_handle &query) const {
                return std::find(vhs.begin(), vhs.end(), query) == vhs.end();
            }

            void insert_constraints(const Vertex_handle &vh, const Vertex_handles &vhs, Building &building) const {

                for (size_t i = 0; i < vhs.size(); ++i)
                    if (vh != vhs[i]) building.cdt.insert_constraint(vh, vhs[i]);
            }

            void collapse_interior_face(const Face_handle &fh, Building &building) const {

                fh->info().color = Color(255, 0, 0);
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H