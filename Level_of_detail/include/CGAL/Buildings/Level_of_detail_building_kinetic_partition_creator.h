#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

// Jean Philippe includes.
#include <CGAL/Buildings/jean_philippe/propagation.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_kinetic_partition_creator {
            
        public:
            using Kernel    = InputKernel;
            using Building  = InputBuilding;
            using Buildings = InputBuildings;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using CDT   = typename Building::CDT;
            using Roofs = typename Building::Roofs;
            
            using Roof = typename Building::Roof;
            using Roof_boundary = typename Roof::Roof_boundary;

            using Building_iterator = typename Buildings::iterator;
            using Building_boundary = typename Building::Boundary;

            using Polygon_boundary = typename Building::Polygon_boundary;
            using Polygons         = typename Building::Polygons;

            using JP_kinetic_propagation = JPTD::Kinetic_Propagation;

            using JP_polyhedron  = JPTD::Partition_Polyhedron;
            using JP_polyhedrons = std::list<JP_polyhedron*>;

            using JP_polyhedrons_iterator = typename JP_polyhedrons::const_iterator;
            using JP_vertex = JPTD::Partition_Vertex;

            using JP_sequence  = std::list<JP_vertex*>;
            using JP_sequences = std::list<JP_sequence>;

            using JP_sequences_iterator = typename JP_sequences::const_iterator;
            using JP_sequence_set = std::set<JP_vertex*>;

            using JP_sequence_iterator     = typename JP_sequence::const_iterator;
            using JP_sequence_set_iterator = typename JP_sequence_set::const_iterator;

            using JP_side  = JPTD::Partition_Side;
            using JP_sides = std::list<JP_side>;

            using JP_sides_iterator = typename JP_sides::const_iterator;
            using JP_facet = JPTD::Partition_Facet;

            using JP_facets = std::list<JP_facet*>;
            using JP_facets_iterator = typename JP_facets::const_iterator;

            using JP_point_3 = JPTD::CGAL_Point_3;
            using JP_conversions = std::map<const JP_vertex*, int>;

            using JP_facet_vertices = std::vector<JP_vertex*>;

            using Polyhedron  = typename Building::Polyhedron;
            using Polyhedrons = typename Building::Polyhedrons;

            using Polyhedron_vertex   = typename Polyhedron::Vertex;
            using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Polyhedron_facet  = typename Polyhedron::Facet;
            using Polyhedron_facets = typename Polyhedron::Facets;

            Level_of_detail_building_kinetic_partition_creator(const CDT &cdt, const FT ground_height) :
            m_cdt(cdt),
            m_ground_height(ground_height)
            { }

            void create_input(Buildings &buildings) const {
                
                if (buildings.size() == 0) return; int count = 0;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit, ++count) {

                    Building &building = bit->second; building.index = count;
					if (building.is_valid) process_building_input(building);
                }
            }

            void create_output(Buildings &buildings) const {
                
                if (buildings.size() == 0) return; int count = 0;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit, ++count) {

                    Building &building = bit->second; building.index = count;
					if (building.is_valid) process_building_output(building);
                }
            }

        private:
            const CDT &m_cdt;
            const FT m_ground_height;
            
            void process_building_input(Building &building) const {
                
                Polygons &polygons = building.polygons;
                set_input(building, polygons);
            }

            void set_input(const Building &building, Polygons &polygons) const {

                polygons.clear();
                const Building_boundary &building_boundary = building.boundaries[0];

                const FT height = FT(2);
                const FT b1 = FT(9) / FT(10);
                const FT b2 = FT(1) - b1;

                // Add boundaries.
                for (size_t i = 0; i < building_boundary.size(); i += 2) {
					
                    const size_t ip = i + 1;
					CGAL_precondition(ip < building_boundary.size());

                    const Point_2 &p = building_boundary[i]->point();
                    const Point_2 &q = building_boundary[ip]->point();

                    const Point_2 a = Point_2(b1 * p.x() + b2 * q.x(), b1 * p.y() + b2 *q.y());
                    const Point_2 b = Point_2(b2 * p.x() + b1 * q.x(), b2 * p.y() + b1 *q.y());

                    const Point_3 p1 = Point_3(a.x(), a.y(), m_ground_height);
                    const Point_3 p2 = Point_3(b.x(), b.y(), m_ground_height);
                    const Point_3 p3 = Point_3(b.x(), b.y(), m_ground_height + height);
                    const Point_3 p4 = Point_3(a.x(), a.y(), m_ground_height + height);

                    Polygon_boundary new_polygon(4);
                    new_polygon[0] = p1;
                    new_polygon[1] = p2;
                    new_polygon[2] = p3;
                    new_polygon[3] = p4;

                    polygons.push_back(new_polygon);
				}

                // Add roofs.
                const Roofs &roofs = building.roofs;
                for (size_t i = 0; i < roofs.size(); ++i) {
					const Roof_boundary &roof_boundary = roofs[i].boundary;

					if (roof_boundary.size() == 0) continue;
					if (!roofs[i].is_valid)        continue;

                    Polygon_boundary new_polygon(roof_boundary.size());
					for (size_t j = 0; j < roof_boundary.size(); ++j) {
							
                        const Point_3 &p = roof_boundary[j];
						new_polygon[j] = p;
					}
                    polygons.push_back(new_polygon);
				}
            }

            void process_building_output(Building &building) const {
                
                const Polygons &polygons = building.polygons;
                JP_kinetic_propagation kinetic(polygons);

	            if (!kinetic.data()) return;
                kinetic.run();

                set_output(kinetic, building);
                kinetic.delete_kinetic_data_structure();
            }

            void set_output(const JP_kinetic_propagation &kinetic, Building &building) const {

                building.polyhedrons.clear();
                for (JP_polyhedrons_iterator p_it = kinetic.partition->polyhedrons_begin(); p_it != kinetic.partition->polyhedrons_end(); ++p_it)
                    add_polyhedron(*p_it, building);

                add_polyhedron_facets(kinetic, building);
            }

            void add_polyhedron(const JP_polyhedron *input_polyhedron, Building &building) const {
                
                Polyhedrons &polyhedrons = building.polyhedrons;
                Polyhedron output_polyhedron;
                
                JP_sequences sequences_per_side; JP_conversions conversions;
                get_polyhedron_vertices(input_polyhedron, output_polyhedron.vertices, sequences_per_side, conversions);
                get_polyhedron_facets(sequences_per_side, conversions, output_polyhedron.facets);

                polyhedrons.push_back(output_polyhedron);
            }

            void get_polyhedron_vertices(const JP_polyhedron *polyhedron, Polyhedron_vertices &vertices, JP_sequences &sequences_per_side, JP_conversions &conversions) const {

                vertices.clear();
                sequences_per_side.clear();
                conversions.clear();

                JP_sequence_set vertices_used;
                for (JP_sides_iterator s_it = polyhedron->facets_begin(); s_it != polyhedron->facets_end(); ++s_it) {
                    const JP_facet *facet = s_it->first;

                    JP_sequence facet_vertices;
                    facet->get_circular_sequence_of_vertices(facet_vertices, !s_it->second);

                    for (JP_sequence_iterator v_it = facet_vertices.begin(); v_it != facet_vertices.end(); ++v_it) vertices_used.insert(*v_it);
                    sequences_per_side.push_back(facet_vertices);
                }

                for (JP_sequence_set_iterator v_it = vertices_used.begin(); v_it != vertices_used.end(); ++v_it) {
                    
                    const JP_vertex  *v = *v_it;
                    const JP_point_3 &p =  v->M;

                    const FT x = static_cast<FT>(CGAL::to_double(p.x()));
                    const FT y = static_cast<FT>(CGAL::to_double(p.y()));
                    const FT z = static_cast<FT>(CGAL::to_double(p.z()));
                    
                    const Polyhedron_vertex vertex = Polyhedron_vertex(x, y, z);

                    vertices.push_back(vertex);
                    conversions[v] = int(vertices.size()) - 1;
                }
            }

            void get_polyhedron_facets(const JP_sequences &sequences_per_side, const JP_conversions &conversions, Polyhedron_facets &facets) const {
                facets.clear();

                for (JP_sequences_iterator s_it = sequences_per_side.begin(); s_it != sequences_per_side.end(); ++s_it) {
                    const JP_sequence &sequence = *s_it;

                    Polyhedron_facet facet;
                    for (JP_sequence_iterator v_it = sequence.begin(); v_it != sequence.end(); ++v_it) facet.push_back(conversions.at(*v_it));
                    facets.push_back(facet);
                }
            }

            void add_polyhedron_facets(const JP_kinetic_propagation &kinetic, Building &building) const {

                Polyhedrons &polyhedrons = building.polyhedron_facets;
                Polyhedron polyhedron;
                
                get_polyhedron_vertices_and_facets(kinetic, polyhedron.vertices, polyhedron.facets);
                polyhedrons.push_back(polyhedron);
            }

            void get_polyhedron_vertices_and_facets(const JP_kinetic_propagation &kinetic, Polyhedron_vertices &vertices, Polyhedron_facets &facets) const {

                JP_facet_vertices v;
                kinetic.partition->octree->get_all_vertices_sorted_by_identifier(v);
                
                CGAL_precondition(v.size() > 0);
                for (size_t i = 0; i < v.size(); ++i) {
                    const JP_point_3 &p = v[i]->M;

                    const FT x = static_cast<FT>(CGAL::to_double(p.x()));
                    const FT y = static_cast<FT>(CGAL::to_double(p.y()));
                    const FT z = static_cast<FT>(CGAL::to_double(p.z()));
                 
                    const Polyhedron_vertex vertex = Polyhedron_vertex(x, y, z);
                    vertices.push_back(vertex);
                }

                const auto &partiton_facets = kinetic.partition->facets;
                for (size_t i = 0; i < partiton_facets.size(); ++i) {
                    for (JP_facets_iterator f_it = partiton_facets[i].begin(); f_it != partiton_facets[i].end(); ++f_it) {
                        
                        const JP_facet* f = *f_it;
                        if (f->p < 6) continue;

                        JP_sequence facet_vertices;
                        f->get_circular_sequence_of_vertices(facet_vertices, true);

                        Polyhedron_facet facet;
                        for (JP_sequence_iterator v_it = facet_vertices.begin(); v_it != facet_vertices.end(); ++v_it) facet.push_back((*v_it)->id);
                        facets.push_back(facet);
                    }
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H