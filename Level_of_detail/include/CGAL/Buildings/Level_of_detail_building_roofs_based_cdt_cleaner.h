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
#include <CGAL/Level_of_detail_enum.h>
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
            using Triangle_2 = typename Kernel::Triangle_2;

            using CDT                = typename Building::CDT;
            using Buildings_iterator = typename Buildings::iterator;

            using Roof          = typename Building::Roof;
            using Roofs         = typename Building::Roofs;
            using Face_boundary = typename Roof::Roof_boundary;

            using Edge                = typename CDT::Edge;
            using Face_handle         = typename CDT::Face_handle;
            using Faces_iterator      = typename CDT::Finite_faces_iterator;
            using CDT_edges_iterator  = typename CDT::Finite_edges_iterator;
            using List_edges          = std::list<Edge>;
            using List_edges_iterator = typename List_edges::iterator;
            
            using Wrong_faces          = std::map<Face_handle, bool>;
            using Wrong_faces_iterator = typename Wrong_faces::const_iterator;

            using Label        = int;
            using Labels       = std::vector<Label>;
            using Sorted_faces = std::map<Face_handle, Labels>;
            using Label_pair   = std::pair<Face_handle, Labels>;
            using Comparator   = std::function<bool(Label_pair, Label_pair)>;
            using Label_set    = std::set<Label_pair, Comparator>;

            using Index   = typename Building::Index;
			using Indices = typename Building::Indices;
            using Shapes  = typename Building::Shapes;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using Plane_associater    = CGAL::LOD::Level_of_detail_building_partition_vote_based_plane_associater<Kernel, Input, Building>;
            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;
            
            using Locate_type = typename CDT::Locate_type;

            Level_of_detail_building_roofs_based_cdt_cleaner(const Input &input, const FT ground_height, Buildings &buildings) : 
            m_input(input),
            m_ground_height(ground_height),
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(1000000)),
            m_thin_face_max_size(FT(1) / FT(2)),
            m_max_percentage(FT(97)) 
            { }

            void clean() {

                m_exterior_faces.clear();
                m_thin_boundary_faces.clear();

                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    remove_exterior_faces(building, true);
                    remove_thin_boundary_faces(building, true);
                    
                    remove_thin_interior_faces(building);
                    remove_faces_with_the_smallest_contribution(building);

                    remove_thin_boundary_faces(building, false);
                    remove_exterior_faces(building, false);

                    update_roofs(building);
                }
            }

        private:
            const Input &m_input;
            const FT m_ground_height;

            Buildings &m_buildings;
            Roof_face_validator m_roof_face_validator;

            const FT m_tolerance;
            const FT m_thin_face_max_size;
            const FT m_max_percentage;

            Wrong_faces m_exterior_faces;
            Wrong_faces m_thin_boundary_faces;

            bool is_exterior_face_1(const Building &building, const Face_handle &fh) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                Face_boundary face_boundary(3);
                face_boundary[0] = Point_3(p1.x(), p1.y(), FT(0));
                face_boundary[1] = Point_3(p2.x(), p2.y(), FT(0));
                face_boundary[2] = Point_3(p3.x(), p3.y(), FT(0));

                return !m_roof_face_validator.is_valid_roof_face(building, face_boundary, true);
            }

            bool is_exterior_face(const Face_handle &fh) const {
                return m_exterior_faces.find(fh) != m_exterior_faces.end();
            }

            bool is_thin_face(const Face_handle &fh, const FT face_size_threshold) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                return triangle.area() < face_size_threshold;
            }

            bool is_thin_boundary_face(const Face_handle &fh) const {
                return m_thin_boundary_faces.find(fh) != m_thin_boundary_faces.end();
            }

            bool is_infinite_face(const CDT &cdt, const Face_handle &fh) const {
                return cdt.is_infinite(fh);
            }

            bool is_boundary_face(const Building &building, const Face_handle &fh) const {

                const CDT &cdt = building.cdt;
                for (size_t i = 0; i < 3; ++i) {

                    const Face_handle &fhn = fh->neighbor(i);
                    if (is_infinite_face(cdt, fhn) || is_exterior_face(fhn)) return true;
                }
                return false;
            }

            bool is_valid_face(const CDT &cdt, const Face_handle &fh) const {
                return !is_infinite_face(cdt, fh) && !is_exterior_face(fh) && !is_thin_boundary_face(fh);
            }

            bool is_exterior_face_2(const Building &building, const Face_handle &fh) const {
                return is_infinite_face(building.cdt, fh) || is_exterior_face_1(building, fh) || (is_thin_face(fh, m_tolerance) && is_boundary_face(building, fh));
            }

            bool can_be_removed(const Face_handle &fh) const {
                return !is_exterior_face(fh) && !is_thin_boundary_face(fh);
            }

            void remove_exterior_faces(Building &building, const bool virtual_deletion) {
                
                CDT &cdt = building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    
                    const Face_handle &fh = static_cast<Face_handle>(fit);
                    if (is_exterior_face_1(building, fh)) {
                        
                        if (virtual_deletion) m_exterior_faces[fh] = true;
                        else cdt.delete_face(fh);
                    }
                }
            }

            void remove_thin_boundary_faces(Building &building, const bool virtual_deletion) {

                CDT &cdt = building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    
                    const Face_handle &fh = static_cast<Face_handle>(fit);
                    if (is_thin_face(fh, m_tolerance) && is_boundary_face(building, fh)) {
                        
                        if (virtual_deletion) m_thin_boundary_faces[fh] = true;
                        else cdt.delete_face(fh);
                    }
                }
            }

            void remove_thin_interior_faces(Building &building) {

                bool state = true;
                CDT &cdt = building.cdt;

                while (state) {

                    state = false;
                    for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

                        const Face_handle &fh = static_cast<Face_handle>(fit);
                        if (can_be_removed(fh) && is_thin_face(fh, m_thin_face_max_size)) {
                            
                            collapse_thin_face(fh, cdt);
                            state = true;
                            break;
                        }
                    }
                }
            }

            void collapse_thin_face(const Face_handle &fh, CDT &cdt) const {

                Edge edge;
                find_smallest_edge(fh, edge);
                cdt.tds().join_vertices(edge);
            }

            void find_smallest_edge(const Face_handle &fh, Edge &edge) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                const FT dist1 = squared_distance_2(p1, p2);
                const FT dist2 = squared_distance_2(p2, p3);
                const FT dist3 = squared_distance_2(p3, p1);

                if (dist1 < dist2 && dist1 < dist3) {
                    edge = std::make_pair(fh, 2); return;
                }
                if (dist2 < dist1 && dist2 < dist3) {
                    edge = std::make_pair(fh, 0); return;
                }
                edge = std::make_pair(fh, 1);
            }

            void remove_faces_with_the_smallest_contribution(Building &building) {

                bool state = true;
                FT current_percentage = FT(100);

                size_t total_size = 0;
                for (size_t i = 0; i < building.shapes.size(); ++i)
                    total_size += building.shapes[i].size();

                while (state) {

                    CDT &cdt = building.cdt;
                    state = false;

                    Sorted_faces sorted_faces;
                    create_sorted_faces(building, sorted_faces);

                    Label_set label_set;
                    sort_faces(sorted_faces, label_set);

                    for (Label_pair label_pair : label_set) {
                        const size_t label_pair_size = compute_label_pair_size(label_pair);
                        
                        const FT percentage     = (static_cast<FT>(label_pair_size) / static_cast<FT>(total_size)) * FT(100);
                        const FT new_percentage = current_percentage - percentage;

                        if (new_percentage >= m_max_percentage) {

                            const Face_handle &fh = label_pair.first;
                            if (!is_exterior_face_2(building, fh)) {
                                current_percentage = new_percentage;

                                std::cout << percentage << ", " << new_percentage << std::endl;

                                const bool status = collapse_boundary_face(fh, building);
                                if (!status) collapse_interior_face(fh, cdt);

                                state = true;
                                break;
                            }
                        } else break;
                    }
                }
            }

            bool collapse_boundary_face(const Face_handle &fh, Building &building) const {

                for (size_t i = 0; i < 3; ++i) {
                    const Face_handle &fhn = fh->neighbor(i);

                    if (is_exterior_face_2(building, fhn)) {

                        const Edge edge = std::make_pair(fh, i);
                        building.cdt.tds().join_vertices(edge);
                        return true;
                    }
                }
                return false;
            }

            void collapse_interior_face(const Face_handle &fh, CDT &cdt) const {
                
                const Edge edge = std::make_pair(fh, 2);
                cdt.tds().join_vertices(edge);

                // collapse_thin_face(fh, cdt);
            }

            void create_sorted_faces(const Building &building, Sorted_faces &sorted_faces) const {

                sorted_faces.clear();

                const CDT &cdt       = building.cdt;
                const Shapes &shapes = building.shapes;

                const size_t num_roof_shapes = shapes.size();
                assert(num_roof_shapes != 0);

                Locate_type locate_type;
				int locate_stub_index = -1;

                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    Labels labels(num_roof_shapes, 0);
                    sorted_faces[fh] = labels;
                }

                for (size_t i = 0; i < num_roof_shapes; ++i) {
                    const Indices &indices = shapes[i];

                    for (size_t j = 0; j < indices.size(); ++j) {
                        const Point_3 &p = m_input.point(indices[j]);

                        const Face_handle fh = cdt.locate(Point_2(p.x(), p.y()), locate_type, locate_stub_index);
					    if (locate_type == CDT::FACE || locate_type == CDT::EDGE || locate_type == CDT::VERTEX)
                            if (!is_exterior_face_2(building, fh)) 
                                sorted_faces.at(fh)[i]++;
                    }
                }
            }

            void sort_faces(const Sorted_faces &sorted_faces, Label_set &label_set) const {

                Comparator comparator = [](const Label_pair &a, const Label_pair &b) {
				    
                    size_t size_a = 0;
                    for (size_t i = 0; i < a.second.size(); ++i)
                        size_a += a.second[i];

                    size_t size_b = 0;
                    for (size_t i = 0; i < b.second.size(); ++i)
                        size_b += b.second[i];
                    
                    return size_a < size_b;
			    };

                label_set = Label_set(sorted_faces.begin(), sorted_faces.end(), comparator);
            }

            size_t compute_total_size(const Label_set &label_set) const {

                size_t total_size = 0;
                for (Label_pair label_pair : label_set)
                    total_size += compute_label_pair_size(label_pair);

                return total_size;
            }

            size_t compute_label_pair_size(const Label_pair &label_pair) const {

                size_t label_pair_size = 0;
                for (size_t i = 0; i < label_pair.second.size(); ++i)
                    label_pair_size += label_pair.second[i];

                return label_pair_size;
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