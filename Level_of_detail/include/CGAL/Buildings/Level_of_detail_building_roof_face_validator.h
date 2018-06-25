#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H

// STL includes.
#include <vector>
#include <cassert>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding>
		class Level_of_detail_building_roof_face_validator {
            
        public:
            typedef KernelTraits  Kernel;
            typedef InputBuilding Building;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_2 = typename Kernel::Triangle_2;
            
            using Boundary = std::vector<Point_3>;

            using CDT            = typename Building::CDT;
            using Face_handle    = typename CDT::Face_handle;
            using Faces_iterator = typename CDT::Finite_faces_iterator;

            using Partition_traits = CGAL::Partition_traits_2<Kernel>;
            using Polygon          = typename Partition_traits::Polygon_2;
            using Polygons         = std::vector<Polygon>;

            Level_of_detail_building_roof_face_validator() :
            m_ghost_tolerance(FT(1) / FT(1000000))
            { }

            bool is_valid_roof_face(const Building &building, const Boundary &boundary, const bool use_barycentre_query_point) const {

                Point_2 query;
                if (use_barycentre_query_point) find_query_point_using_barycentre(boundary, query);
                else find_query_point_using_partition(boundary, query);

                const auto &faces = building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                        
                    const Point_2 &p1 = faces[i]->vertex(0)->point();
                    const Point_2 &p2 = faces[i]->vertex(1)->point();
                    const Point_2 &p3 = faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (triangle.has_on_bounded_side(query) || triangle.has_on_boundary(query)) return true;
                }
                return false;
            }

            bool is_ghost_face(const Face_handle &fh) const {
                
                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                return triangle.area() < m_ghost_tolerance;
            }

            bool is_outside_face(const Building &building, const Face_handle &fh) const {

                const Point_2 &p1 = fh->vertex(0)->point();
                const Point_2 &p2 = fh->vertex(1)->point();
                const Point_2 &p3 = fh->vertex(2)->point();

                Boundary face_boundary(3);
                face_boundary[0] = Point_3(p1.x(), p1.y(), FT(0));
                face_boundary[1] = Point_3(p2.x(), p2.y(), FT(0));
                face_boundary[2] = Point_3(p3.x(), p3.y(), FT(0));

                return !is_valid_roof_face(building, face_boundary, true);
            }

            void mark_wrong_faces(Building &building) const {
                
                CDT &cdt = building.cdt;
                for (Faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
                    const Face_handle &fh = static_cast<Face_handle>(fit);

                    if (fit->info().is_checked) continue;
                    if (is_ghost_face(fit)) fit->info().is_valid = false;
                    if (is_outside_face(building, fit)) fit->info().is_valid = false;
                }
            }

        private:
            const FT m_ghost_tolerance;
            
            void find_query_point_using_barycentre(const Boundary &boundary, Point_2 &query) const {

                FT x = FT(0), y = FT(0);
                for (size_t i = 0; i < boundary.size(); ++i) {
                    const Point_3 &p = boundary[i];
                    
                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(boundary.size());
                y /= static_cast<FT>(boundary.size());

                query = Point_2(x, y);
            }

            void find_query_point_using_partition(const Boundary &boundary, Point_2 &query) const {
                
                const auto &points = boundary;
                assert(points.size() > 0);

                Polygon polygon;
                for (size_t i = 0; i < points.size(); ++i)
                    polygon.push_back(Point_2(points[i].x(), points[i].y()));

                if (orientation_2(polygon.vertices_begin(), polygon.vertices_end(), Partition_traits()) != CGAL::COUNTERCLOCKWISE) 
                    polygon.reverse_orientation();

                Polygons result;
                CGAL::approx_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(), std::back_inserter(result));

                FT x = FT(0), y = FT(0);
                for (auto vit = result[0].vertices_begin(); vit != result[0].vertices_end(); ++vit) {

                    x += (*vit).x();
                    y += (*vit).y();
                }

                const FT n = static_cast<FT>(std::distance(result[0].vertices_begin(), result[0].vertices_end()));

                x /= n;
                y /= n;

                query = Point_2(x, y);
            }
        };
        
    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H