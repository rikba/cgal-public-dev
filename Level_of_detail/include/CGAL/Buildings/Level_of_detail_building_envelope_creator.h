#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_CREATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/Cartesian.h>
#include <CGAL/envelope_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_building_roof_face_validator.h>
#include <CGAL/Buildings/Associaters/Level_of_detail_building_envelope_plane_associater.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputBuildings, class InputBuilding>
		class Level_of_detail_building_envelope_creator {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT         = typename Kernel::FT;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_3 = typename Kernel::Triangle_3;

            using Building_iterator = typename Buildings::iterator;

            using Roof              = typename Building::Roof;
            using Data              = typename Building::Data;
            using Envelope_input    = typename Building::Data_triangles;
            using Envelope_element  = typename Building::Data_triangle;
            using Associated_planes = typename Roof::Associated_planes;

            using Exact_type       = CGAL::Exact_rational;
            using Exact_kernel     = CGAL::Cartesian<Exact_type>;
            using Env_traits       = CGAL::Env_triangle_traits_3<Exact_kernel>;
            using Exact_point      = typename Exact_kernel::Point_3;
            using Exact_triangle   = typename Env_traits::Surface_3;
            using Data_traits      = CGAL::Env_surface_data_traits_3<Env_traits, Data>;
            using Data_triangle    = typename Data_traits::Surface_3;
            using Data_triangles   = std::vector<Data_triangle>;
            using Envelope_diagram = CGAL::Envelope_diagram_2<Data_traits>;

            using Face_iterator       = typename Envelope_diagram::Face_const_iterator;
            using Halfedge_circulator = typename Envelope_diagram::Ccb_halfedge_const_circulator;
            
            using Boundary = std::vector<Point_3>;

            using Plane_associater    = CGAL::LOD::Level_of_detail_building_envelope_plane_associater<Kernel, Building, Envelope_diagram>;
            using Roof_face_validator = CGAL::LOD::Level_of_detail_building_roof_face_validator<Kernel, Building>;

            Level_of_detail_building_envelope_creator(const Input &, const FT ground_height) :
            m_ground_height(ground_height) 
            { }

            void create(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    Building &building = bit->second;
					if (building.is_valid) process_building(building);
                }
            }

            void set_min_face_width(const FT) { }

        private:
            const FT m_ground_height;
            
            Plane_associater    m_plane_associater;
            Roof_face_validator m_roof_face_validator;

            void process_building(Building &building) const {
                
                const Envelope_input &envelope_input = building.envelope_input;
                if (envelope_input.size() < 2) {
                    
                    building.is_valid = false;
                    return;
                }

                Data_triangles triangles;
                set_envelope_input(envelope_input, triangles);
                
                Envelope_diagram diagram;
                compute_envelope(triangles, diagram);

                update_roofs(diagram, triangles, building);
            }

            void set_envelope_input(const Envelope_input &envelope_input, Data_triangles &triangles) const {
                for (size_t i = 0; i < envelope_input.size(); ++i) {
                    
                    const Envelope_element &element = envelope_input[i];
                    add_element(element, triangles);
                }
            }

            void add_element(const Envelope_element &element, Data_triangles &triangles) const {
                
                const Triangle_3 &tri = element.first;
                const Data &data      = element.second;

                const Point_3 &p1 = tri.vertex(0);
                const Point_3 &p2 = tri.vertex(1);
                const Point_3 &p3 = tri.vertex(2);

                const Exact_point ep1 = Exact_point(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p1.z()));
                const Exact_point ep2 = Exact_point(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), CGAL::to_double(p2.z()));
                const Exact_point ep3 = Exact_point(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()), CGAL::to_double(p3.z()));

                const Exact_triangle etri = Exact_triangle(ep1, ep2, ep3);
                const Data_triangle  dtri = Data_triangle(etri, data);
                
                triangles.push_back(dtri);
            }

            void compute_envelope(const Data_triangles &triangles, Envelope_diagram &diagram) const {
                CGAL::lower_envelope_3(triangles.begin(), triangles.end(), diagram);
            }

            void update_roofs(const Envelope_diagram &diagram, const Data_triangles &triangles, Building &building) const {
                
                Roof roof;
                building.clear_roofs();
                
                Boundary &boundary = roof.boundary;
                Associated_planes &associated_planes = roof.associated_planes;

                for (Face_iterator fit = diagram.faces_begin(); fit != diagram.faces_end(); ++fit) {
                    
                    boundary.clear();
                    associated_planes.clear();
                    
                    if (fit->is_unbounded()) continue;

                    Halfedge_circulator ccb = fit->outer_ccb();
                    do {

                        const FT x = static_cast<FT>(CGAL::to_double(ccb->target()->point().x()));
                        const FT y = static_cast<FT>(CGAL::to_double(ccb->target()->point().y()));
                        const FT z = building.height + m_ground_height;

                        boundary.push_back(Point_3(x, y, z));
                        ++ccb;

                    } while (ccb != fit->outer_ccb());

                    if (m_roof_face_validator.is_valid_roof_face(building, boundary)) {
                        
                        m_plane_associater.find_associated_planes(fit, associated_planes);
                        building.roofs.push_back(roof);
                    }
                }
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_CREATOR_H