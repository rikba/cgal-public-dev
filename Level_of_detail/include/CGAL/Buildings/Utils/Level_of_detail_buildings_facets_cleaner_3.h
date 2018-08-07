#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

// New CGAL includes.
#include <CGAL/Region_growing/Level_of_detail_planar_region_growing.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputHDS, class InputBuilding, class FacetHandle>
		class Local_builder : public Modifier_base<InputHDS> {
		
		public:
			using Kernel   	   = InputKernel;
    		using HDS      	   = InputHDS;
            using Building     = InputBuilding;
            using Facet_handle = FacetHandle;

			using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

			using Builder = CGAL::Polyhedron_incremental_builder_3<HDS>;
            using Facets  = std::vector< std::vector<Facet_handle> >;

			Local_builder(const Clean_facets &clean_facets) : 
			m_clean_facets(clean_facets), 
    		m_index_counter(0) { }

			void operator()(HDS &hds) {
				
                m_facets.clear();
                m_facets.resize(1);

				Builder builder(hds, false);

				const size_t expected_num_vertices = estimate_number_of_vertices();
				const size_t expected_num_facets   = estimate_number_of_facets();

				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_facets);
				add_clean_facets(builder);
				builder.end_surface();
			}

            void get_facets(Facets &facets) const {
                facets = m_facets;
            }

		private:
			const Clean_facets &m_clean_facets;
            Facets m_facets;
			size_t m_index_counter;

			inline size_t estimate_number_of_vertices() {
				return m_clean_facets.size() * 4;
			}

			inline size_t estimate_number_of_facets() {
				return m_clean_facets.size();
			}

			void add_clean_facets(Builder &builder) {

				for (size_t i = 0; i < m_clean_facets.size(); ++i)
					add_clean_facet(m_clean_facets[i], builder);
			}

			void add_clean_facet(const Clean_facet &clean_facet, Builder &builder) {

				if (clean_facet.size() == 0) return;

				for (size_t i = 0; i < clean_facet.size(); ++i) 
					builder.add_vertex(clean_facet[i]);

				const Facet_handle fh = builder.begin_facet();
				for (size_t i = 0; i < clean_facet.size(); ++i) builder.add_vertex_to_facet(m_index_counter++);
				builder.end_facet();

                m_facets[0].push_back(fh);
			}
        };

		template<class InputKernel, class InputBuilding, class InputBuildings>
		class Level_of_detail_buildings_facets_cleaner_3 {
            
        public:
            typedef InputKernel    Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;

            using Buildings_iterator = typename Buildings::iterator;

            using Polyhedron   = typename Building::Polyhedron;
			using Polyhedrons  = typename Building::Polyhedrons;
			using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

			using Polyhedron_facet    = typename Polyhedron::Facet;
			using Polyhedron_facets   = typename Polyhedron::Facets;
			using Polyhedron_vertices = typename Polyhedron::Vertices;

            using Mesh = CGAL::Polyhedron_3<Kernel>;
            
            using Mesh_facet_handle   = typename Mesh::Facet_const_handle;
            using Mesh_facets         = std::vector< std::vector<Mesh_facet_handle> >;
            using HDS                 = typename Mesh::HalfedgeDS;
            using Halfedge_handle     = typename Mesh::Halfedge_const_handle;
            using Facet_vertex_handle = typename Mesh::Facet::Vertex_const_handle;

            using Local_builder         = Local_builder<Kernel, HDS, Building, Mesh_facet_handle>;
            using Planar_region_growing = CGAL::LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Mesh_facets>;

            using Building_region  = std::vector<Mesh_facet_handle>;
			using Building_regions = std::vector<Building_region>;
			using Regions 		   = std::vector<Building_regions>;

            Level_of_detail_buildings_facets_cleaner_3(Buildings &buildings) :
            m_buildings(buildings),
            m_tolerance(FT(1) / FT(10000))
            { }

            void create_clean_facets() {
                
                if (m_buildings.size() == 0) 
                    return;

				for (Buildings_iterator b_it = m_buildings.begin(); b_it != m_buildings.end(); ++b_it) {
                    Building &building = b_it->second;

					if (building.is_valid) 
                        process_building(building);
                }
            }

        private:
            Buildings &m_buildings;
            const FT m_tolerance;
            
            void process_building(Building &building) const {
            
                Clean_facets &clean_facets = building.clean_facets;
                clean_facets.clear();
                
                building.is_clean = true;
                const Polyhedrons &polyhedrons = building.polyhedrons;
                
				for (size_t i = 0; i < polyhedrons.size(); ++i)
					process_polyhedron(i, polyhedrons, clean_facets);

                merge_clean_facets(clean_facets);
            }

            void process_polyhedron(const size_t polyhedron_index, const Polyhedrons &polyhedrons, Clean_facets &clean_facets) const {

                const Polyhedron &polyhedron = polyhedrons[polyhedron_index];
                Clean_facet clean_facet;

                if (!polyhedron.is_valid) return;

				const Polyhedron_facets   &facets   = polyhedron.facets;
				const Polyhedron_vertices &vertices = polyhedron.vertices;

                for (size_t i = 0; i < facets.size(); ++i) {
                    const Polyhedron_facet &facet = facets[i];

                    if (!facet.is_valid) continue;

                    clean_facet.clear();
                    clean_facet.resize(facet.indices.size());

                    for (size_t j = 0; j < facet.indices.size(); ++j)
                        clean_facet[j] = vertices[facet.indices[j]];

                    if (!is_interior_facet(clean_facet, polyhedron_index, polyhedrons)) clean_facets.push_back(clean_facet);
                }
            }

            bool is_interior_facet(const Clean_facet &clean_facet, const size_t polyhedron_index, const Polyhedrons &polyhedrons) const {

                for (size_t i = 0; i < polyhedrons.size(); ++i) {
                    const Polyhedron &polyhedron = polyhedrons[i];
                    
                    if (!polyhedron.is_valid || i == polyhedron_index) continue;

                    const Polyhedron_facets   &facets   = polyhedron.facets;
				    const Polyhedron_vertices &vertices = polyhedron.vertices;

                    for (size_t j = 0; j < facets.size(); ++j) {
                        if (!facets[j].is_valid) continue;

                        if (are_equal_facets(clean_facet, facets[j], vertices)) 
                            return true;
                    }
                }
                return false;
            }

            bool are_equal_facets(const Clean_facet &clean_facet, const Polyhedron_facet &facet, const Polyhedron_vertices &vertices) const {

                if (clean_facet.size() != facet.indices.size()) return false;

                size_t count = 0;
                for (size_t i = 0; i < clean_facet.size(); ++i) {
                    for (size_t j = 0; j < facet.indices.size(); ++j) {
                        
                        if (are_equal_points(clean_facet[i], vertices[facet.indices[j]])) {
                            
                            ++count;
                            break;
                        }
                    }
                }
                return count == clean_facet.size();
            }

            bool are_equal_points(const Point_3 &p, const Point_3 &q) const {
                return CGAL::abs(p.x() - q.x()) < m_tolerance && CGAL::abs(p.y() - q.y()) < m_tolerance && CGAL::abs(p.z() - q.z()) < m_tolerance;
            }

            void merge_clean_facets(Clean_facets &clean_facets) const {

                Mesh mesh;

                Local_builder builder(clean_facets);
                mesh.delegate(builder);

                Mesh_facets facets;
                builder.get_facets(facets);

                Planar_region_growing planar_region_growing(facets);

				Regions regions;
				planar_region_growing.find_regions(regions);

                set_merged_facets(regions[0], clean_facets);
            }

            void set_merged_facets(const Building_regions &regions, Clean_facets &clean_facets) const {
                clean_facets.clear();

                Clean_facet clean_facet;
                for (size_t i = 0; i < regions.size(); ++i) {

                    create_clean_facet(regions[i], clean_facet);
                    clean_facets.push_back(clean_facet);
                }
            }

            void create_clean_facet(const Building_region &region, Clean_facet &clean_facet) const {
                clean_facet.clear();
                
                if (region.size() == 1) {

                    const Mesh_facet_handle &fh = region[0];
                    add_clean_facet(fh, clean_facet);
                }

                Halfedge_handle curr, end; size_t region_index = 1000000;
                find_first_halfedge_handle(region, curr, region_index);

                // clean_facet.push_back(curr->vertex()->point());
                // end = curr;

                // curr = curr->next();
                // while (curr != end) {

                //     find_next_halfedge_handle(region, region_index, curr, clean_facet);
                // };
            }

            void add_clean_facet(const Mesh_facet_handle &fh, Clean_facet &clean_facet) const {

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                clean_facet.push_back(vh->point());

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    clean_facet.push_back(vh->point());

                    curr = curr->next();
                };
            }

            void find_first_halfedge_handle(const Building_region &region, Halfedge_handle &he, size_t &region_index) const {

                for (size_t i = 0; i < region.size(); ++i) {
                    region_index = i;

                    const Mesh_facet_handle &fh = region[region_index];
                    
                    he = fh->halfedge();
                    if (is_valid(he, region_index, region)) return;

                    Halfedge_handle end = he;
                    he = he->next();

                    while (he != end) {

                        if (is_valid(he, region_index, region)) return;
                        he = he->next();
                    }
                }
            }

            bool is_valid(const Halfedge_handle &he, const size_t region_index, const Building_region &region) const {

                const Facet_vertex_handle &vh = he->vertex();
                const Point_3 &p = vh->point();

                Halfedge_handle stub;
                for (size_t i = 0; i < region.size(); ++i) {
                    if (i == region_index) continue;
                    
                    const Mesh_facet_handle &fh = region[i];
                    if (contains(p, fh, stub)) return false;
                }
                return true;
            }

            bool contains(const Point_3 &p, const Mesh_facet_handle &fh, Halfedge_handle &next) const {

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                if (are_equal_points(p, vh->point())) {
                 
                    next = curr;
                    return true;
                }

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    if (are_equal_points(p, vh->point())) {
                     
                        next = curr;
                        return true;
                    }

                    curr = curr->next();
                };
                return false;
            }

            void find_next_halfedge_handle(const Building_region &region, size_t &region_index, Halfedge_handle &curr, Clean_facet &clean_facet) const {

                const Point_3 &p = curr->vertex()->point();
                if (is_valid(curr, region_index, region)) {

                    clean_facet.push_back(p);
                    curr = curr->next();
                    
                    return;
                }

                Halfedge_handle next;
                for (size_t i = 0; i < region.size(); ++i) {
                    if (i == region_index) continue;
                    
                    const Mesh_facet_handle &fh = region[i];
                    if (contains(p, fh, next)) {

                        clean_facet.push_back(p);
                        curr = next->next();

                        region_index = i;
                        return;
                    }
                }

                curr = curr->next();
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_FACETS_CLEANER_3_H