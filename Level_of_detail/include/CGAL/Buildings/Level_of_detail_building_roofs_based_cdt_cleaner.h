#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H

// STL includes.
#include <map>
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Buildings/Level_of_detail_building_roof_face_validator.h>
#include <CGAL/Buildings/Associaters/Level_of_detail_building_partition_vote_based_plane_associater.h>

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
            using Face_circulator = typename CDT::Face_circulator;

            using Roof    = typename Building::Roof;
            using Roofs   = typename Building::Roofs;
            using Indices = typename Building::Indices;

            using Contributions             = typename Building::Contributions;
            using Face_contributions        = typename Building::Face_contributions;
            using Comparator                = typename Building::Comparator;
            using Contribution_pair         = typename Building::Contribution_pair;
            using Sorted_face_contributions = typename Building::Sorted_face_contributions;

            using Locate_type = typename CDT::Locate_type;
            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using Plane_associater    = CGAL::LOD::Level_of_detail_building_partition_vote_based_plane_associater<Kernel, Input, Building>;
            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;

            using Color = CGAL::Color;
            using Log   = CGAL::LOD::Mylog;
            
            using Vertex_handles = std::vector<Vertex_handle>;
            using Ids            = std::vector< std::pair<int, bool> >;
            using Edges          = std::vector<Edge>;

            Level_of_detail_building_roofs_based_cdt_cleaner(const Input &input, const FT ground_height, Buildings &buildings) : 
            m_input(input),
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_thin_face_max_size(FT(1) / FT(2)),
            m_max_percentage(FT(90)),
            m_max_main_iters(50),
            m_max_circulator_iters(30),
            m_angle_tolerance(FT(1))
            { }

            void clean() {

                size_t count = 0; // std::cout << std::endl;
                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit, ++count) {
                    Building &building = (*bit).second;
                    
                    // std::cout << FT(count) / FT(m_buildings.size()) * FT(100) << " %" << std::endl;
                    clean_cdt(building);

                    m_roof_face_validator.mark_wrong_faces(building);
                    update_roofs(building);
                }
                // std::cout << FT(count) / FT(m_buildings.size()) * FT(100) << " %" << std::endl;
            }

        private:
            const Input &m_input;
            const FT m_ground_height;

            Buildings &m_buildings;

            const FT m_thin_face_max_size;
            const FT m_max_percentage;

            Roof_face_validator m_roof_face_validator;

            const size_t m_max_main_iters;
            const size_t m_max_circulator_iters;

            const FT m_angle_tolerance;

            void clean_cdt(Building &building) const {

                building.total_contributions_size = building.interior_indices.size();
                bool collapse_detected = false; size_t iters = 0;

                do {
                    building.current_percentage = FT(100);

                    // Default indices.
                    CDT &cdt = building.cdt; int count = 0; ++iters;
                    for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
                        fit->info().index = count++;

                    // Contributions.
                    create_face_contributions(building);
                    sort_face_contributions(building);
                    finilize_face_contributions(building);

                    // Handle faces.
                    for (auto sc_it = building.sorted_face_contributions.begin(); sc_it != building.sorted_face_contributions.end(); ++ sc_it) {
                        Face_handle fh = static_cast<Face_handle>(sc_it->first);

                        if (!m_roof_face_validator.is_valid_face(building, fh)) continue;
                        collapse_detected = handle_face(fh, building);

                        if (collapse_detected) {
                         
                            m_roof_face_validator.mark_wrong_faces(building);
                            break;
                        }
                    }

                    // Handle faces.
                    // for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    //     Face_handle fh = static_cast<Face_handle>(fit);

                    //     if (!m_roof_face_validator.is_valid_face(building, fh)) continue;
                    //     collapse_detected = handle_face(fh, building);

                    //     if (collapse_detected) {
                         
                    //         m_roof_face_validator.mark_wrong_faces(building);
                    //         break;
                    //     }
                    // }

                } while (collapse_detected && iters < m_max_main_iters);

                if (iters >= m_max_main_iters) {
                 
                    // std::cout << "Exceeded number of main iterations!" << std::endl;
                    return;
                }
            }

            void create_face_contributions(Building &building) const {
                
                Face_contributions &face_contributions = building.face_contributions;
                face_contributions.clear();

                const CDT &cdt         = building.cdt;
                const Indices &indices = building.interior_indices;

                const size_t num_roof_shapes = 1;
                CGAL_precondition(num_roof_shapes != 0);

                Locate_type locate_type;
				int locate_stub_index = -1;

                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    Contributions contributions(num_roof_shapes, 0);
                    face_contributions[fh] = contributions;
                }

                for (size_t i = 0; i < indices.size(); ++i) {
                    const Point_3 &p = m_input.point(indices[i]);

                    const Face_handle fh = cdt.locate(Point_2(p.x(), p.y()), locate_type, locate_stub_index);
					if (locate_type == CDT::FACE || locate_type == CDT::EDGE || locate_type == CDT::VERTEX)
                        face_contributions.at(fh)[0]++;
                }
            }

            void sort_face_contributions(Building &building) const {
                
                Face_contributions        &face_contributions        = building.face_contributions;
                Sorted_face_contributions &sorted_face_contributions = building.sorted_face_contributions;

                Comparator comparator = [](const Contribution_pair &a, const Contribution_pair &b) {
				    
                    size_t size_a = 0;
                    for (size_t i = 0; i < a.second.size(); ++i)
                        size_a += a.second[i];

                    size_t size_b = 0;
                    for (size_t i = 0; i < b.second.size(); ++i)
                        size_b += b.second[i];
                    
                    return size_a < size_b;
			    };

                sorted_face_contributions.clear();
                sorted_face_contributions.resize(face_contributions.size());

                size_t i = 0;
                for (auto fc_it = face_contributions.begin(); fc_it != face_contributions.end(); ++fc_it, ++i)
                    sorted_face_contributions[i] = *fc_it;

                face_contributions.clear();
                std::sort(sorted_face_contributions.begin(), sorted_face_contributions.end(), comparator);
            }

            void finilize_face_contributions(Building &building) const {

                FT &current_percentage = building.current_percentage;
                const size_t &total_contributions_size = building.total_contributions_size;

                Sorted_face_contributions &sorted_face_contributions = building.sorted_face_contributions;                
                for (auto sc_it = sorted_face_contributions.begin(); sc_it != sorted_face_contributions.end(); ++sc_it) {
                    
                    Face_handle fh = sc_it->first;
                    if (!m_roof_face_validator.is_valid_face(building, fh)) continue;

                    const size_t face_contributions_size = compute_face_contribution_size(*sc_it);
                    const FT percentage = (static_cast<FT>(face_contributions_size) / static_cast<FT>(total_contributions_size)) * FT(100);

                    const FT new_percentage = current_percentage - percentage;
                    // std::cout << "percent = " << current_percentage << "," << new_percentage << ", " << face_contributions_size << std::endl;

                    if (new_percentage >= m_max_percentage) {
                        
                        current_percentage = new_percentage;
                        continue;

                    } else {

                        sorted_face_contributions.erase(sc_it, sorted_face_contributions.end());
                        break;
                    }
                }
            }

            size_t compute_face_contribution_size(const Contribution_pair &contribution_pair) const {

                size_t contribution_pair_size = 0;
                for (size_t i = 0; i < contribution_pair.second.size(); ++i)
                    contribution_pair_size += contribution_pair.second[i];

                return contribution_pair_size;
            }

            bool handle_face(const Face_handle &fh, Building &building) const {
                
                /* if (should_be_collapsed(building, fh)) */ return collapse_face(fh, building);
                // return false;
            }

            bool should_be_collapsed(const Building &building, const Face_handle &fh) const {
                
                if (is_insignificant(building, fh)) return true;
                if (is_thin(fh)) return true;

                return false;
            }

            inline bool is_thin(const Face_handle &fh) const {
                return m_roof_face_validator.is_ghost_face(fh, m_thin_face_max_size);
            }

            bool is_insignificant(const Building &building, const Face_handle &fh) const {
                const Sorted_face_contributions &sorted_face_contributions = building.sorted_face_contributions;

                for (auto sc_it = sorted_face_contributions.begin(); sc_it != sorted_face_contributions.end(); ++sc_it) {
                    if (sc_it->first == fh) {
                        return true;
                    }
                }
                return false;
            }

            bool collapse_face(const Face_handle &fh, Building &building) const {

                Ids ids;
                check_boundary_face_edges(fh, building, ids);

                CGAL_precondition(ids.size() >= 0 && ids.size() <= 3);
                if (ids.size() == 1) {
                    
                    if (!ids[0].second) return false;
                    return collapse_boundary_face(fh, ids[0].first, building);
                }
                else if (ids.size() == 2) return false;
                else if (ids.size() == 3) return false;
                
                return collapse_interior_face(fh, building);
            }

            void check_boundary_face_edges(const Face_handle &fh, const Building &building, Ids &ids) const {
                
                ids.clear();
                for (int i = 0; i < 3; ++i) {
                    Face_handle fhn = fh->neighbor(i);
                    
                    if (!m_roof_face_validator.is_valid_face(building, fhn)) {
                        
                        const bool is_valid_edge = check_boundary_edge_validity(building, fh, i);
                        ids.push_back(std::make_pair(i, is_valid_edge));
                    }
                }
            }

            bool check_boundary_edge_validity(const Building &building, const Face_handle &fh, const int vertex_index) const {

                const Edge edge = Edge(fh, vertex_index);
                
                Vertex_handle vh1, vh2;
                find_edge_vertices(edge, vh1, vh2);

                if (is_corner_vertex(building, vh1, fh)) return false;
                if (is_corner_vertex(building, vh2, fh)) return false;

                return true;
            }

            bool check_interior_edge_validity(const Building &building, const Face_handle &fh, const int vertex_index) const {

                const Edge edge = Edge(fh, vertex_index);
                
                Vertex_handle vh1, vh2;
                find_edge_vertices(edge, vh1, vh2);

                if (is_boundary_vertex(building, vh1, fh)) return false;
                if (is_boundary_vertex(building, vh2, fh)) return false;

                return true;
            }

            FT compute_vertex_angle(const Building &building, const Vertex_handle &vh, const Face_handle &fh) const {

                const CDT &cdt = building.cdt;

                Face_circulator face_circulator = cdt.incident_faces(vh, fh);
                Face_circulator end = face_circulator;

                if (face_circulator.is_empty()) return -FT(1);
                size_t iters = 0;

                FT total_angle = FT(0);
                do {

                    Face_handle face = static_cast<Face_handle>(face_circulator);
                    if (!m_roof_face_validator.is_valid_face(building, face)) {

                        ++face_circulator;
                        ++iters;    
                        
                        if (face_circulator == end || iters >= m_max_circulator_iters) break;
                        else continue;
                    }

                    const FT angle = compute_angle(vh, face);
                    total_angle += angle;
                    
                    ++face_circulator;
                    ++iters;

                } while (face_circulator != end && iters < m_max_circulator_iters);
                
                if (iters >= m_max_circulator_iters) {
                    
                    std::cout << "Exceeded number of face circulator iterations!" << std::endl;
                    return -FT(1);
                }

                return total_angle;
            }

            FT compute_angle(const Vertex_handle &vh, const Face_handle &fh) const {

                const int curr = fh->index(vh);
                const int prev = (curr + 2) % 3;
                const int next = (curr + 1) % 3;

                const Point_2 &a = fh->vertex(prev)->point();
                const Point_2 &b = fh->vertex(curr)->point();
                const Point_2 &c = fh->vertex(next)->point();

                const Point_2 ab = Point_2(b.x() - a.x(), b.y() - a.y());
                const Point_2 cb = Point_2(b.x() - c.x(), b.y() - c.y());

                const FT dot   = ab.x() * cb.x() + ab.y() * cb.y();
                const FT cross = ab.x() * cb.y() - ab.y() * cb.x();

                const FT alpha = static_cast<FT>(atan2(CGAL::to_double(cross), CGAL::to_double(dot)));

                return CGAL::abs(static_cast<FT>(floor(CGAL::to_double(alpha * FT(180) / CGAL_PI + FT(1) / FT(2)))));
            }

            bool is_corner_vertex(const Building &building, const Vertex_handle &vh, const Face_handle &fh) const {

                const FT total_angle = compute_vertex_angle(building, vh, fh);
                if (total_angle < FT(0)) return true;

                // std::cout << "angle = " << total_angle << std::endl;

                if (CGAL::abs(total_angle - FT(180)) <= m_angle_tolerance) 
                    return false;

                if (CGAL::abs(total_angle - FT(360)) <= m_angle_tolerance) 
                    return false;
                
                return true;
            }

            bool is_boundary_vertex(const Building &building, const Vertex_handle &vh, const Face_handle &fh) const {

                const FT total_angle = compute_vertex_angle(building, vh, fh);
                if (total_angle < FT(0)) return true;

                // std::cout << "angle = " << total_angle << std::endl;

                if (CGAL::abs(total_angle - FT(180)) <= m_angle_tolerance) 
                    return true;

                if (CGAL::abs(total_angle - FT(360)) <= m_angle_tolerance) 
                    return false;
                
                return true;
            }

            inline bool collapse_boundary_face(const Face_handle &fh, const int vertex_index, Building &building) const {
                
                // std::cout << "boundary ";
                const Edge boundary_edge = std::make_pair(fh, vertex_index);
                return collapse_edge(boundary_edge, building);
            }

            bool collapse_edge(const Edge &edge, Building &building) const {
                CGAL_precondition(edge.second >= 0 && edge.second < 3);

                Vertex_handle vh1, vh2;
                find_edge_vertices(edge, vh1, vh2);

                Point_2 new_point;
                find_new_edge_point(vh1, vh2, new_point);

                const Face_handle &fh = edge.first;

                Vertex_handles vhs1;
                bool success = update_incident_constraints(vh1, fh, vhs1, building);

                if (!success) return false;

                Vertex_handles vhs2;
                success = update_incident_constraints(vh2, fh, vhs2, building);

                if (!success) return false;

                Vertex_handles vhs;
                merge_constraints(vhs1, vhs2, vh1, vh2, vhs);

                building.cdt.remove(vh1);
                building.cdt.remove(vh2);

                const Vertex_handle vh = building.cdt.insert(new_point);
                insert_constraints(vh, vhs, building);

                // std::cout << "edge collapse finished" << std::endl;
                return true;
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

            bool update_incident_constraints(const Vertex_handle &query, const Face_handle &fh, Vertex_handles &vhs, Building &building) const {
                
                vhs.clear();
                CDT &cdt = building.cdt;

                Edge_circulator edge_circulator = cdt.incident_edges(query, fh);
                Edge_circulator end = edge_circulator;

                if (edge_circulator.is_empty()) return false;
                size_t iters = 0;

                do {
                    ++edge_circulator;
                    ++iters;
                } while (edge_circulator != end && iters < m_max_circulator_iters);

                if (iters >= m_max_circulator_iters) {
                    
                    std::cout << "Exceeded number of edge circulator iterations!" << std::endl;
                    return false;
                }

                Edges edges;
                edge_circulator = end; iters = 0;

                do {
                    const Edge &edge = *edge_circulator;
                    if (cdt.is_constrained(edge) && !cdt.is_infinite(edge)) {

                        Vertex_handle vh1, vh2;
                        find_edge_vertices(edge, vh1, vh2);

                        CGAL_precondition(vh1 == query || vh2 == query);
                        if (vh1 != query) vhs.push_back(vh1);
                        if (vh2 != query) vhs.push_back(vh2);

                        edges.push_back(edge);
                    }
                    
                    ++edge_circulator;
                    ++iters;

                } while (edge_circulator != end && iters < m_max_circulator_iters);

                for (size_t i = 0; i < edges.size(); ++i)
                    cdt.remove_constrained_edge(edges[i].first, edges[i].second);

                return true;
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

            bool collapse_interior_face(const Face_handle &fh, Building &building) const {

                // return false;
                Edge smallest_edge = std::make_pair(fh, -1);

                find_smallest_edge(fh, smallest_edge);
                const int index1 = smallest_edge.second;

                CGAL_precondition(smallest_edge.second != -1);
                bool is_valid_edge = check_interior_edge_validity(building, smallest_edge.first, smallest_edge.second);

                if (is_valid_edge) {
                    
                    // std::cout << "interior ";
                    return collapse_edge(smallest_edge, building);

                } else {
                    
                    find_smallest_edge(fh, smallest_edge);
                    const int index2 = smallest_edge.second;

                    CGAL_precondition(smallest_edge.second != -1);
                    is_valid_edge = check_interior_edge_validity(building, smallest_edge.first, smallest_edge.second);

                    if (is_valid_edge) {
                        
                        return collapse_edge(smallest_edge, building);
                        // std::cout << "interior ";

                    } else {

                        return false;

                        /*
                        CGAL_precondition(smallest_edge.second != -1);
                        if ( (index1 == 0 && index2 == 1) || (index1 == 1 && index2 == 0) ) {
                         
                            smallest_edge = std::make_pair(fh, 2);
                            return collapse_edge(smallest_edge, building);
                        }

                        if ( (index1 == 1 && index2 == 2) || (index1 == 2 && index2 == 1) ) {
                         
                            smallest_edge = std::make_pair(fh, 0);
                            return collapse_edge(smallest_edge, building);
                        }

                        smallest_edge = std::make_pair(fh, 1);
                        return collapse_edge(smallest_edge, building); */
                    }
                }
                return false;
            }

            void find_smallest_edge(const Face_handle &fh, Edge &edge) const {

                const Point_2 &p0 = fh->vertex(0)->point();
                const Point_2 &p1 = fh->vertex(1)->point();
                const Point_2 &p2 = fh->vertex(2)->point();

                const FT dist0 = squared_distance_2(p0, p1);
                const FT dist1 = squared_distance_2(p1, p2);
                const FT dist2 = squared_distance_2(p2, p0);

                if (edge.second == -1) {

                    if (dist0 < dist1 && dist0 < dist2) {
                        edge = std::make_pair(fh, 2); return;
                    }
                    if (dist1 < dist2 && dist1 < dist0) {
                        edge = std::make_pair(fh, 0); return;
                    }
                    edge = std::make_pair(fh, 1); return;
                }

                if (edge.second == 0) {

                    if (dist0 < dist2) {
                        edge = std::make_pair(fh, 2); return;
                    }
                    edge = std::make_pair(fh, 1); return;
                }

                if (edge.second == 1) {

                    if (dist0 < dist1) {
                        edge = std::make_pair(fh, 2); return;
                    }
                    edge = std::make_pair(fh, 0); return;
                }

                if (edge.second == 2) {

                    if (dist1 < dist2) {
                        edge = std::make_pair(fh, 0); return;
                    }
                    edge = std::make_pair(fh, 1); return;
                }
            }

            void update_roofs(Building &building) const {

                add_roof_faces(building);
                add_associated_planes(building);
            }

            void add_roof_faces(Building &building) const {
                
                const CDT &cdt = building.cdt;
                const FT  z    = building.height + m_ground_height;
                
                Roofs &roofs = building.roofs;
                roofs.clear();

                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    if (!fit->info().is_valid) continue;
                    
                    Roof roof;
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    const Point_2 &p1 = fh->vertex(0)->point();
                    const Point_2 &p2 = fh->vertex(1)->point();
                    const Point_2 &p3 = fh->vertex(2)->point();

                    roof.boundary.push_back(Point_3(p1.x(), p1.y(), z));
                    roof.boundary.push_back(Point_3(p2.x(), p2.y(), z));
                    roof.boundary.push_back(Point_3(p3.x(), p3.y(), z));

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

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H