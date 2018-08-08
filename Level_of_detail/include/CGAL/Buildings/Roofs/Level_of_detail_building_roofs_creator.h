#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H

// STL includes.
#include <set>
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

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
            using Point_3 = typename Kernel::Point_3;
            
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

            using Unique_facets = std::set<Facet>;

            Level_of_detail_building_roofs_creator(const Input &input, const FT ground_height, Buildings &buildings) :
            m_input(input),
            m_buildings(buildings),
            m_ground_height(ground_height)
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

            void process_building(Building &building) const {

                Roofs &roofs = building.roofs;
    
                const Polyhedrons &polyhedrons = building.polyhedrons;
                const Shapes &shapes           = building.shapes;

                for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    create_roofs(indices, polyhedrons, roofs);
                }
            }

            void create_roofs(const Indices &indices, const Polyhedrons &polyhedrons, Roofs &roofs) const {

                Facets facets;
                find_closest_facets(indices, polyhedrons, facets);

                remove_duplicates(facets);
                pick_best_coverage_facets(indices, facets);

                create_roofs_from_facets(facets, roofs);
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

            }

            void create_facet(const Polyhedron_facet &poly_facet, const Polyhedron_vertices &poly_vertices, Facet &facet) const {
                
                facet.clear();
                facet.resize(poly_facet.indices.size());

                for (size_t i = 0; i < poly_facet.indices.size(); ++i)
                    facet[i] = poly_vertices[poly_facet.indices[i]];
            }

            void remove_duplicates(Facets &facets) const {

                Unique_facets unique_facets;
                create_unique_facets(facets, unique_facets);
                set_facets_from_unique_facets(unique_facets, facets);
            }

            void create_unique_facets(const Facets &facets, Unique_facets &unique_facets) const {

            }

            void set_facets_from_unique_facets(const Unique_facets &unique_facets, Facets &facets) const {
                
                facets.clear();
                facets.resize(unique_facets.size());

                size_t i = 0;
                for (auto uf_it = unique_facets.begin(); uf_it != unique_facets.end(); ++uf_it, ++i)
                    facets[i] = *uf_it;
            }

            void pick_best_coverage_facets(const Indices &indices, Facets &facets) const {

                sort_facets_by_size(facets);
                pick_best_coverage(indices, facets);
            }

            void sort_facets_by_size(Facets &facets) const {

            }

            void pick_best_coverage(const Indices &indices, Facets &facets) const {

            }

            void create_roofs_from_facets(const Facets &facets, Roofs &roofs) const {
                
                Roof roof; roofs.clear();
                for (size_t i = 0; i < facets.size(); ++i) {
                    
                    roof.boundary = facets[i];
                    roofs.push_back(roof);
                }
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOFS_CREATOR_H