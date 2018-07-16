#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_FACE_VALIDATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_FACE_VALIDATOR_H

// CGAL includes.
#include <CGAL/utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding>
		class Level_of_detail_building_cdt_face_validator {
            
        public:
            typedef InputKernel   Kernel;
            typedef InputBuilding Building;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            using CDT            = typename Building::CDT;
            using Face_handle    = typename CDT::Face_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;

            Level_of_detail_building_cdt_face_validator(Building &building) :
            m_building(building),
            m_ghost_tolerance(FT(1) / FT(1000000))
            { }

            bool is_ghost_face(const Face_handle &fh, const FT ghost_tolerance) const {
                
                Triangle_2 triangle;
                construct_triangle(fh, triangle);

                return triangle.area() < ghost_tolerance;
            }

            bool is_ghost_boundary_face(const Face_handle &fh) const {
                return is_ghost_face(fh, m_ghost_tolerance) && is_boundary_face(fh);
            }

            bool is_boundary_face(const Face_handle &fh) const {
                
                for (size_t i = 0; i < 3; ++i) {
                    const Face_handle &fhn = fh->neighbor(i);
                    
                    if (m_building.cdt.is_infinite(fhn) || is_exterior_face(fhn)) 
                        return true;
                }
                return false;
            }

            bool is_exterior_face(const Face_handle &fh) const {

                Point_2 query;
                construct_query_point(fh, query);

                const auto &faces = m_building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                    
                    Triangle_2 triangle;
                    construct_triangle(faces[i], triangle);

                    if (triangle.has_on_bounded_side(query) || triangle.has_on_boundary(query)) return false;
                }
                return true;
            }

            bool is_interior_face(const Face_handle &fh) const {
                return !is_exterior_face(fh);
            }

            void validate(Face_handle &fh) const {

                if (m_building.cdt.is_infinite(fh)) {
                    fh->info().is_valid = false; return;
                }

                if (is_ghost_boundary_face(fh)) {
                    fh->info().is_valid = false; return;
                }

                if (is_exterior_face(fh)) {
                    fh->info().is_valid = false; return;
                }

                fh->info().is_valid = true; return;
            }

            void validate() const {
                
                CDT &cdt = m_building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    
                    Face_handle fh = static_cast<Face_handle>(fit);
                    validate(fh);
                }
            }

        private:
            Building &m_building;
            const FT m_ghost_tolerance;
            
            void construct_triangle(const Face_handle &fh, Triangle_2 &triangle) const {
                
                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                triangle = Triangle_2(p1, p2, p3);
            }
            
            void construct_query_point(const Face_handle &fh, Point_2 &query) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                const FT x = (p1.x() + p2.x() + p3.x()) / FT(3);
                const FT y = (p1.y() + p2.y() + p3.y()) / FT(3);

                query = Point_2(x, y);
            }
        };
        
    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_CDT_FACE_VALIDATOR_H