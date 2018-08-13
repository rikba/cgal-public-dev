#ifndef CGAL_LEVEL_OF_DETAIL_LOCAL_MESH_BUILDER_H
#define CGAL_LEVEL_OF_DETAIL_LOCAL_MESH_BUILDER_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputHDS, class InputBuilding, class FacetHandle>
		class Local_mesh_builder : public Modifier_base<InputHDS> {
		
		public:
			using Kernel   	   = InputKernel;
    		using HDS      	   = InputHDS;
            using Building     = InputBuilding;
            using Facet_handle = FacetHandle;

			using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

			using Builder = CGAL::Polyhedron_incremental_builder_3<HDS>;
            using Facets  = std::vector< std::vector<Facet_handle> >;

			Local_mesh_builder(const Clean_facets &clean_facets) : 
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

				if (clean_facet.first.size() == 0) return;

				for (size_t i = 0; i < clean_facet.first.size(); ++i) 
					builder.add_vertex(clean_facet.first[i]);

				const Facet_handle fh = builder.begin_facet();
				for (size_t i = 0; i < clean_facet.first.size(); ++i) builder.add_vertex_to_facet(m_index_counter++);
				builder.end_facet();

                m_facets[0].push_back(fh);
			}
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_LOCAL_MESH_BUILDER_H