#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H

// STL includes.
#include <vector>
#include <utility>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputContainer, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roofs_creator {
            
        public:
            using Kernel    = InputKernel;
            using Input     = InputContainer;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;
            
            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Index   = typename Building::Index;
            using Indices = typename Building::Indices;
            using Shapes  = typename Building::Shapes;
            using Roof    = typename Building::Roof;
            using Roofs   = typename Building::Roofs;

            using Polyhedron_facet    = typename Polyhedron::Facet;
            using Polyhedron_facets   = typename Polyhedron::Facets;
            using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Facet  = std::vector<Point_3>;
            using Facets = std::vector<Facet>;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            using Polygon       = CGAL::Polygon_2<Kernel>;
            using Pair          = std::pair<size_t, size_t>;
            using Polygon_pair  = std::pair<Polygon, Pair>;
            using Polygon_pairs = std::vector<Polygon_pair>;

            using Coordinates = std::vector<FT>;

            using Mean_value = CGAL::Barycentric_coordinates::Mean_value_2<Kernel>;
            using Mean_value_coordinates = CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel>;

            using Hull = std::vector<Point_2>;

            Level_of_detail_building_roofs_creator(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height),
            m_distance_threshold(FT(1) / FT(2)),
            m_tolerance(FT(1) / FT(100000)),
            m_polygon_area_tolerance(FT(1) / FT(10000))
            { }

            void create_roofs() {

                if (m_buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid && building.shapes.size() != 0 && building.jp_polygons.size() != 0)
                        process_building(building);
                }
            }

        private:
            const Input &m_input;
            Buildings   &m_buildings;

            const FT m_ground_height;
            const FT m_distance_threshold;

            const FT m_tolerance;
            const FT m_polygon_area_tolerance;

            void process_building(Building &building) const {

                Roofs &roofs             = building.roofs;
                Polyhedrons &polyhedrons = building.polyhedrons;

                const Shapes &shapes = building.shapes;
                preprocess(roofs, polyhedrons);

                for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    add_roofs(indices, polyhedrons, roofs);
                }
            }

            void preprocess(Roofs &roofs, Polyhedrons &polyhedrons) const {

                roofs.clear();
            }

            void add_roofs(const Indices &indices, const Polyhedrons &polyhedrons, Roofs &roofs) const {

                Facets facets;
                find_closest_facets(indices, polyhedrons, facets);

                remove_duplicates(facets);
                pick_best_coverage_facets(indices, facets);

                add_roofs_from_facets(facets, roofs);
            }

            void find_closest_facets(const Indices &indices, const Polyhedrons &polyhedrons, Facets &facets) const {

                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    const Polyhedron &polyhedron = polyhedrons[i];

                    process_polyhedron(indices, polyhedron, facets);
                }
            }

            void process_polyhedron(const Indices &indices, const Polyhedron &polyhedron, Facets &facets) const {

                const Polyhedron_facets   &poly_facets   = polyhedron.facets;
                const Polyhedron_vertices &poly_vertices = polyhedron.vertices;
                
                for (size_t i = 0; i < poly_facets.size(); ++i) {
                    const Polyhedron_facet &poly_facet = poly_facets[i];

                    if (is_close_poly_facet(indices, poly_facet, poly_vertices)) {

                        Facet facet;
                        create_facet(poly_facet, poly_vertices, facet);
                        facets.push_back(facet);
                    }
                }
            }

            bool is_close_poly_facet(const Indices &indices, const Polyhedron_facet &poly_facet, const Polyhedron_vertices &poly_vertices) const {

                Plane_3 plane;
                create_plane(poly_facet, poly_vertices, plane);

                const FT average_squared_distance = compute_average_squared_distance_to_plane(indices, plane);

                if (average_squared_distance < FT(0)) return false;
                return average_squared_distance < m_distance_threshold * m_distance_threshold;
            }

            void create_plane(const Polyhedron_facet &poly_facet, const Polyhedron_vertices &poly_vertices, Plane_3 &plane) const {

                CGAL_precondition(poly_facet.indices.size() > 2);

                const Point_3 &p1 = poly_vertices[poly_facet.indices[0]];
                const Point_3 &p2 = poly_vertices[poly_facet.indices[1]];
                const Point_3 &p3 = poly_vertices[poly_facet.indices[2]];

                plane = Plane_3(p1, p2, p3);
            }

            FT compute_average_squared_distance_to_plane(const Indices &indices, const Plane_3 &plane) const {

                if (indices.size() == 0) return -FT(1);

                FT average_squared_distance = FT(0);
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const Point_3 &p = m_input.point(indices[i]);
                    const Point_3  q = plane.projection(p);

                    const FT squared_distance = squared_distance_3(p, q);
                    average_squared_distance += squared_distance;
                }
                
                average_squared_distance /= static_cast<FT>(indices.size());
                return average_squared_distance;
            }

            void create_facet(const Polyhedron_facet &poly_facet, const Polyhedron_vertices &poly_vertices, Facet &facet) const {
                
                facet.clear();
                facet.resize(poly_facet.indices.size());

                for (size_t i = 0; i < poly_facet.indices.size(); ++i)
                    facet[i] = poly_vertices[poly_facet.indices[i]];
            }

            void remove_duplicates(Facets &facets) const {

                Facets new_facets;
                for (size_t i = 0; i < facets.size(); ++i) {
                    const Facet &facet = facets[i];

                    if (!contains(facet, new_facets)) 
                        new_facets.push_back(facet);
                }
                
                facets.clear();
                facets = new_facets;
            }

            bool contains(const Facet &query_facet, const Facets &facets) const {

                for (size_t i = 0; i < facets.size(); ++i)
                    if (are_aqual_facets(query_facet, facets[i])) return true;
                return false;
            }

            bool are_aqual_facets(const Facet &f1, const Facet &f2) const {
                
                if (f1.size() != f2.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < f1.size(); ++i) {
                    for (size_t j = 0; j < f2.size(); ++j) {

                        if (CGAL::abs(f1[i].x() - f2[j].x()) < m_tolerance && CGAL::abs(f1[i].y() - f2[j].y()) < m_tolerance && CGAL::abs(f1[i].z() - f2[j].z()) < m_tolerance) {
                            ++count; break;
                        }
                    }
                }

                return count == f1.size();
            }

            void pick_best_coverage_facets(const Indices &indices, Facets &facets) const {

                Polygon_pairs polygon_pairs;
                
                Hull hull;
                compute_convex_hull(indices, hull);
                
                create_polygon_pairs(facets, polygon_pairs);
                compute_polygon_inliners(indices, polygon_pairs);
                
                sort_polygon_pairs(polygon_pairs);
                pick_best_coverage(polygon_pairs, hull, facets);
            }

            void compute_convex_hull(const Indices &indices, Hull &hull) const {
                hull.clear();

				Hull points(indices.size());
				for (size_t i = 0; i < indices.size(); ++i) 
                    points[i] = Point_2(m_input.point(indices[i]).x(), m_input.point(indices[i]).y());

				CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));
            }

            void create_polygon_pairs(const Facets &facets, Polygon_pairs &polygon_pairs) const {
                
                polygon_pairs.clear();
                polygon_pairs.resize(facets.size());

                for (size_t i = 0; i < facets.size(); ++i)
                    create_polygon_pair(facets[i], i, polygon_pairs[i]);
            }

            void create_polygon_pair(const Facet &facet, const size_t facet_index, Polygon_pair &polygon_pair) const {

                Polygon polygon;
                create_polygon(facet, polygon);

                polygon_pair = std::make_pair(polygon, std::make_pair(facet_index, 0));
            }

            void create_polygon(const Facet &facet, Polygon &polygon) const {

                polygon.clear();
                for (size_t i = 0; i < facet.size(); ++i) {
                    const Point_3 &p = facet[i];

                    const Point_2 q = Point_2(p.x(), p.y());
                    polygon.push_back(q);
                }

                if (polygon.is_clockwise_oriented()) polygon.reverse_orientation();
            }

            void compute_polygon_inliners(const Indices &indices, Polygon_pairs &polygon_pairs) const {

                for (size_t i = 0; i < indices.size(); ++i) {
                    const Point_3 &p = m_input.point(indices[i]);

                    const Point_2 q = Point_2(p.x(), p.y());
                    add_inliers(q, polygon_pairs);
                }
            }

            void add_inliers(const Point_2 &query, Polygon_pairs &polygon_pairs) const {

                for (size_t i = 0; i < polygon_pairs.size(); ++i) {
                    Polygon_pair &polygon_pair = polygon_pairs[i];
                    
                    const Polygon &polygon = polygon_pair.first;
                    if (!is_valid_polygon(polygon)) continue;

                    if (belongs(query, polygon.vertices_begin(), polygon.vertices_end())) {

                        Pair &pair = polygon_pair.second;
                        pair.second++;

                        return;
                    }
                }
            }

            bool is_valid_polygon(const Polygon &polygon) const {

                if (polygon.size() < 3)                        return false;
                if (polygon.area() < m_polygon_area_tolerance) return false;
                if (!polygon.is_simple())                      return false;

                return true;
            }

            template<class Vertices_iterator>
            bool belongs(const Point_2 &query, const Vertices_iterator &vertices_begin, const Vertices_iterator &vertices_end) const {

                Mean_value_coordinates mean_value_coordinates(vertices_begin, vertices_end);

                Coordinates coordinates;
                mean_value_coordinates(query, std::back_inserter(coordinates));

                for (size_t i = 0; i < coordinates.size(); ++i)
                    if (coordinates[i] < FT(0)) return false;
                return true;
            }

            void sort_polygon_pairs(Polygon_pairs &polygon_pairs) const {

                std::sort(polygon_pairs.begin(), polygon_pairs.end(), [](const Polygon_pair &p1, const Polygon_pair &p2) -> bool {
                    
                    const Pair &pair1 = p1.second;
                    const Pair &pair2 = p2.second;
                    
                    return pair1.second > pair2.second;
                });
            }

            void pick_best_coverage(const Polygon_pairs &polygon_pairs, const Hull &hull, Facets &facets) const {

                Facets new_facets;
                FT num_points = FT(0);

                for (size_t i = 0; i < polygon_pairs.size(); ++i) {
                    
                    const Polygon_pair &polygon_pair = polygon_pairs[i];
                    const Pair &pair = polygon_pair.second;

                    const size_t facet_index = pair.first;

                    if (pair.second != 0 && should_be_included(polygon_pair.first, hull))
                        new_facets.push_back(facets[facet_index]);
                }

                facets.clear();
                facets = new_facets;
            }

            bool should_be_included(const Polygon &polygon, const Hull &hull) const {

                Point_2 b1;
                compute_barycentre(polygon.vertices_begin(), polygon.vertices_end(), b1);

                // Point_2 b2;
                // compute_barycentre(hull.begin(), hull.end(), b2);

                return belongs(b1, hull.begin(), hull.end()); // && belongs(b2, polygon.vertices_begin(), polygon.vertices_end());
            }

            template<class Vertices_iterator>
            void compute_barycentre(const Vertices_iterator &vertices_begin, const Vertices_iterator &vertices_end, Point_2 &b) const {

                FT x = FT(0), y = FT(0), size = FT(0);
                for (auto v_it = vertices_begin; v_it != vertices_end; ++v_it, size += FT(1)) {
                    const Point_2 &p = *v_it;

                    x += p.x();
                    y += p.y();
                }

                x /= size;
                y /= size;

                b = Point_2(x, y);
            }

            void add_roofs_from_facets(const Facets &facets, Roofs &roofs) const {
                
                Roof roof;
                for (size_t i = 0; i < facets.size(); ++i) {
                    
                    roof.boundary = facets[i];
                    roofs.push_back(roof);
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H