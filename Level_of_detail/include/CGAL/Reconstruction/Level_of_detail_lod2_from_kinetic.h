#ifndef CGAL_LEVEL_OF_DETAIL_LOD2_FROM_KINETIC_H
#define CGAL_LEVEL_OF_DETAIL_LOD2_FROM_KINETIC_H

// STL includes.
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputHDS, class InputCDT, class InputBuilding, class InputBuildings, class FacetHandle>
		class LOD2_builder_from_kinetic : public Modifier_base<InputHDS> {
		
		public:
			using Kernel   			 = InputKernel;
    		using HDS      			 = InputHDS;
			using CDT      			 = InputCDT;
    		using Building 		     = InputBuilding;
			using Buildings 		 = InputBuildings;
    		using Color_facet_handle = FacetHandle;
    		using Facet_handle 		 = FacetHandle;

    		using FT 	     = typename Kernel::FT;
			using Point_2    = typename Kernel::Point_2;
			using Point_3    = typename Kernel::Point_3;
			using Triangle_2 = typename Kernel::Triangle_2;

			using Boundary 	  		= typename Building::Boundary;			
			using Floor_faces 		= typename Building::Floor_faces;
			using Floor_face_handle = typename Building::Face_handle;
			using Vertex_handle		= typename Building::Vertex_handle;
			using Polyhedron  		= typename Building::Polyhedron;
			using Polyhedrons 	    = typename Building::Polyhedrons;
			using Clean_facet       = typename Building::Clean_facet;
			using Clean_facets      = typename Building::Clean_facets;

			using Polyhedron_facet    	  = typename Polyhedron::Facet;
			using Polyhedron_facets   	  = typename Polyhedron::Facets;
			using Polyhedron_vertices 	  = typename Polyhedron::Vertices;

			using Color 	   = CGAL::Color;
			using Facet_colors = std::map<Color_facet_handle, Color>;
			using Builder 	   = CGAL::Polyhedron_incremental_builder_3<HDS>;
    		
			using Building_iterator = typename Buildings::const_iterator;
    		using Ground 			= std::vector<Point_3>;

			LOD2_builder_from_kinetic(const CDT &cdt, const Buildings &buildings, const Ground &ground, const FT ground_height, Facet_colors &facet_colors) : 
    		m_cdt(cdt),
			m_buildings(buildings), 
			m_ground(ground),
			m_ground_height(ground_height),
			m_facet_colors(facet_colors),
    		m_index_counter(0) { }

			void operator()(HDS &hds) {
				
				Builder builder(hds, false);

				const size_t expected_num_vertices = estimate_number_of_vertices();
				const size_t expected_num_facets   = estimate_number_of_facets();

				m_index_counter = 0;
				builder.begin_surface(expected_num_vertices, expected_num_facets);

				for (Building_iterator bit = m_buildings.begin(); bit != m_buildings.end(); ++bit) {
					const Building &building = (*bit).second;

					if (!building.is_valid) 
						continue;

					if (building.shapes.size() == 0 || building.jp_polygons.size() == 0) {
						
						add_new_lod1_building(building, builder);
						continue;
					}

					add_new_lod2_building(building, builder);
				}

				add_ground(builder);
				builder.end_surface();
			}

		private:
			const CDT 		&m_cdt;
			const Buildings &m_buildings;
			const Ground    &m_ground;
			
			const FT m_ground_height;
			Facet_colors &m_facet_colors;
			
			size_t m_index_counter;

			inline size_t estimate_number_of_vertices() {
				return m_buildings.size() * 8 + 4;
			}

			inline size_t estimate_number_of_facets() {
				return m_buildings.size() * 6 + 1;
			}

			void add_new_lod1_building(const Building &building, Builder &builder) {

				const FT height   = building.height;
				const Color color = building.color;

				const Floor_faces &floor_faces = building.faces;
				add_horizontal_triangulation(floor_faces, color, height, builder);
	
				const Boundary &boundary = building.boundaries[0];
				add_walls_from_unoriented_boundary(boundary, color, FT(0), height, builder);
			}

			void add_horizontal_triangulation(const Floor_faces &floor_faces, const Color &color, const FT height, Builder &builder) {

				const size_t num_faces = floor_faces.size();
				for (size_t i = 0; i < num_faces; ++i) {

					const Floor_face_handle &fh = floor_faces[i];
					const Triangle_2 &triangle  = m_cdt.triangle(fh);

					add_triangle_face(triangle, color, height, builder);
				}
			}

			void add_triangle_face(const Triangle_2 &triangle, const Color &color, const FT height, Builder &builder) {

				const Point_2 &va = triangle.vertex(0);
				const Point_2 &vb = triangle.vertex(1);
				const Point_2 &vc = triangle.vertex(2);

				const Point_3 a = Point_3(va.x(), va.y(), m_ground_height + height);
				const Point_3 b = Point_3(vb.x(), vb.y(), m_ground_height + height);
				const Point_3 c = Point_3(vc.x(), vc.y(), m_ground_height + height);

				add_triangle(a, b, c, color, builder);
			}

			void add_triangle(const Point_3 &a, const Point_3 &b, const Point_3 &c, const Color &color, Builder &builder) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);

		        const Color_facet_handle cfh = builder.begin_facet();

		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        
		        builder.end_facet();
		        m_facet_colors[cfh] = color;
			}

			void add_walls_from_unoriented_boundary(const Boundary &boundary, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {
				
				const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; i += 2) {
					
					const size_t ip = i + 1;
					assert(ip < num_vertices);

					add_quadrilateral_wall(boundary[i], boundary[ip], color, height_floor, height_roof, builder);
				}
			}

			void add_quadrilateral_wall(const Vertex_handle &vi, const Vertex_handle &vj, const Color &color, const FT height_floor, const FT height_roof, Builder &builder) {

				const Point_2 &va = vi->point();
				const Point_2 &vb = vj->point();

				const Point_3 a = Point_3(va.x(), va.y(), m_ground_height + height_floor);
				const Point_3 b = Point_3(vb.x(), vb.y(), m_ground_height + height_floor);

				const Point_3 c = Point_3(vb.x(), vb.y(), m_ground_height + height_roof);
				const Point_3 d = Point_3(va.x(), va.y(), m_ground_height + height_roof);

				add_quad(a, b, c, d, color, builder);
			}

			void add_quad(const Point_3 &a, const Point_3 &b, const Point_3 &c, const Point_3 &d, const Color &color, Builder &builder) {

		        builder.add_vertex(a);
		        builder.add_vertex(b);
		        builder.add_vertex(c);
		        builder.add_vertex(d);

		        const Color_facet_handle cfh = builder.begin_facet();

		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);
		        builder.add_vertex_to_facet(m_index_counter++);

		        builder.end_facet();
		        m_facet_colors[cfh] = color;
			}

			void add_new_lod2_building(const Building &building, Builder &builder) {

				if (building.is_clean) {
					
					use_clean_version(building, builder);
					return;
				}

				const Color color = building.color;
				const Polyhedrons &polyhedrons = building.polyhedrons;

				for (size_t i = 0; i < polyhedrons.size(); ++i) {
					const Polyhedron &polyhedron = polyhedrons[i];

					add_polyhedron(polyhedron, color, builder);
				}
			}

			void add_polyhedron(const Polyhedron &polyhedron, const Color &color, Builder &builder) {

				if (!polyhedron.is_valid) 
					return;
				
				const Polyhedron_facets   &facets   = polyhedron.facets;
				const Polyhedron_vertices &vertices = polyhedron.vertices;

				for (size_t i = 0; i < facets.size(); ++i)
					add_polyhedron_facet(facets[i], vertices, color, builder);
			}

			void add_polyhedron_facet(const Polyhedron_facet &facet, const Polyhedron_vertices &vertices, const Color &color, Builder &builder) {
				
				size_t count = 0;
				for (size_t i = 0; i < facet.indices.size(); ++i) {
					
					if (!facet.is_valid) {

						++count;
						continue;
					}

					const Point_3 &p = vertices[facet.indices[i]];
					builder.add_vertex(p);
				}

				if (count == facet.indices.size()) return;
				const Color_facet_handle cfh = builder.begin_facet();
				
				for (size_t i = 0; i < facet.indices.size(); ++i) {
					if (!facet.is_valid) continue;

					builder.add_vertex_to_facet(m_index_counter++);
				}
				builder.end_facet();
		        		
				m_facet_colors[cfh] = color;
			}

			void use_clean_version(const Building &building, Builder &builder) {

				const Color 	   &color 		 = building.color;
				const Clean_facets &clean_facets = building.clean_facets;

				for (size_t i = 0; i < clean_facets.size(); ++i)
					add_clean_facet(clean_facets[i], color, builder);
			}

			void add_clean_facet(const Clean_facet &clean_facet, const Color &color, Builder &builder) {

				if (clean_facet.first.size() == 0) return;

				for (size_t i = 0; i < clean_facet.first.size(); ++i) 
					builder.add_vertex(clean_facet.first[i]);

				const Color_facet_handle cfh = builder.begin_facet();
				for (size_t i = 0; i < clean_facet.first.size(); ++i) builder.add_vertex_to_facet(m_index_counter++);
				builder.end_facet();

				m_facet_colors[cfh] = color; // clean_facet.second;
			}

			void add_ground(Builder &builder) {
				assert(!m_ground.empty());

				const size_t num_vertices = m_ground.size();
				assert(num_vertices == 4);

				const Point_3 &a = m_ground[3];
				const Point_3 &b = m_ground[2];
				const Point_3 &c = m_ground[1];
				const Point_3 &d = m_ground[0];

				const Point_3 p1 = Point_3(a.x(), a.y(), a.z() + m_ground_height);
				const Point_3 p2 = Point_3(b.x(), b.y(), b.z() + m_ground_height);
				const Point_3 p3 = Point_3(c.x(), c.y(), c.z() + m_ground_height);
				const Point_3 p4 = Point_3(d.x(), d.y(), d.z() + m_ground_height);

				const Color color(169, 169, 169);
				add_quad(p1, p2, p3, p4, color, builder);
			}
        };

        template<class InputKernel, class InputCDT, class InputBuilding, class InputBuildings, class OutputMesh>
		class Level_of_detail_lod2_from_kinetic {

		public:
			using Kernel    = InputKernel;
			using CDT 		= InputCDT;
			using Building  = InputBuilding;
			using Buildings = InputBuildings;
			using Mesh	    = OutputMesh;

			using FT  = typename Kernel::FT;
            using HDS = typename Mesh::HalfedgeDS;

            using Mesh_facet_handle = typename Mesh::Facet_const_handle;
            using Mesh_builder 		= LOD2_builder_from_kinetic<Kernel, HDS, CDT, Building, Buildings, Mesh_facet_handle>;
			using Mesh_facet_colors = typename Mesh_builder::Facet_colors;
			using Ground 			= typename Mesh_builder::Ground;

            Level_of_detail_lod2_from_kinetic(const CDT &cdt, const Buildings &buildings, const Ground &ground, const FT ground_height, Mesh_facet_colors &mesh_facet_colors) :
            m_builder(cdt, buildings, ground, ground_height, mesh_facet_colors) { }

            inline void reconstruct(Mesh &mesh) {
				mesh.delegate(m_builder);
            }

        private:
			Mesh_builder m_builder;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_LOD2_FROM_KINETIC_H