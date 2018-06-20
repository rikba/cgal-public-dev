#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H

// STL includes.
#include <list>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_building_roof_face_validator.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_based_cdt_cleaner {

        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_2 = typename Kernel::Triangle_2;

            using CDT                = typename Building::CDT;
            using Buildings_iterator = typename Buildings::iterator;

            using Roof          = typename Building::Roof;
            using Face_boundary = typename Roof::Roof_boundary;

            using Edge                = typename CDT::Edge;
            using Face_handle         = typename CDT::Face_handle;
            using Faces_iterator      = typename CDT::Finite_faces_iterator;
            using CDT_edges_iterator  = typename CDT::Finite_edges_iterator;
            using List_edges          = std::list<Edge>;
            using List_edges_iterator = typename List_edges::iterator;

            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;

            Level_of_detail_building_roofs_based_cdt_cleaner(Buildings &buildings) : 
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(1000000))
            { }

            void clean() {

                for (Buildings_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
                    Building &building = (*bit).second;

                    remove_exterior_faces(building);

                    List_edges edges;
                    create_edges(building, edges);
                    update_faces(edges, building);
                }
            }

        private:
            Buildings &m_buildings;
            Roof_face_validator m_roof_face_validator;

            const FT m_tolerance;

            void remove_exterior_faces(Building &building) const {
                
                CDT &cdt = building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    
                    const Face_handle &fh = static_cast<Face_handle>(fit);
                    if (is_exterior_face(building, fh)) cdt.delete_face(fh);
                }
            }

            bool is_exterior_face(const Building &building, const Face_handle &fh) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                Face_boundary face_boundary(3);
                face_boundary[0] = Point_3(p1.x(), p1.y(), FT(0));
                face_boundary[1] = Point_3(p2.x(), p2.y(), FT(0));
                face_boundary[2] = Point_3(p3.x(), p3.y(), FT(0));

                return !m_roof_face_validator.is_valid_roof_face(building, face_boundary, true);
            }

            void create_edges(const Building &building, List_edges &edges) const {
                
                edges.clear();
                const CDT &cdt = building.cdt;

                for (CDT_edges_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
                    
                    const Edge &edge = *eit;
                    if (is_valid_edge(cdt, edge)) edges.push_back(edge);
                }
            }

            bool is_valid_edge(const CDT &cdt, const Edge &edge) const {

                const Face_handle &fh_original = edge.first;
                const Face_handle &fh_opposite = fh_original->neighbor(edge.second);

                const Point_2 &p1 = fh_original->vertex(0)->point();
                const Point_2 &p2 = fh_original->vertex(1)->point();
                const Point_2 &p3 = fh_original->vertex(2)->point();

                const Triangle_2 triangle1 = Triangle_2(p1, p2, p3);

                const Point_2 &q1 = fh_opposite->vertex(0)->point();
                const Point_2 &q2 = fh_opposite->vertex(1)->point();
                const Point_2 &q3 = fh_opposite->vertex(2)->point();

                const Triangle_2 triangle2 = Triangle_2(q1, q2, q3);

                return !cdt.is_infinite(fh_original)  && 
                       !cdt.is_infinite(fh_opposite)  && 
                       triangle1.area() > m_tolerance && 
                       triangle2.area() > m_tolerance;
            }

            void update_faces(List_edges &edges, Building &building) const {
                CDT &cdt = building.cdt;

                for (List_edges_iterator eit = edges.begin(); eit != edges.end(); ++eit) {
                    Edge &edge = *eit;

                    
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_BASED_CDT_CLEANER_H