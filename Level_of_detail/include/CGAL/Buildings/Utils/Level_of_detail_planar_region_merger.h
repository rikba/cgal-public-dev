#ifndef CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H
#define CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <cassert>
#include <stdlib.h>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polyhedron_3.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Buildings/Utils/Level_of_detail_local_mesh_builder.h>
#include <CGAL/Region_growing/Level_of_detail_planar_region_growing.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputBuilding>
		class Level_of_detail_planar_region_merger {

		public:
			using Kernel   = InputKernel;
            using Building = InputBuilding;

			using FT       = typename Kernel::FT;
            using Point_2  = typename Kernel::Point_2;
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;

			using CDT          = typename Building::CDT;
            using Clean_facet  = typename Building::Clean_facet;
			using Clean_facets = typename Building::Clean_facets;

            using Mesh = CGAL::Polyhedron_3<Kernel>;

            using HDS                 = typename Mesh::HalfedgeDS;
            using Mesh_facet_handle   = typename Mesh::Facet_const_handle;
            using Halfedge_handle     = typename Mesh::Halfedge_const_handle;
            using Facet_vertex_handle = typename Mesh::Facet::Vertex_const_handle;

            using Input_region  = std::vector<Mesh_facet_handle>;
            using Input_regions = std::vector<Input_region>;
            using Mesh_facets   = Input_regions;

            using Local_builder         = CGAL::LOD::Local_mesh_builder<Kernel, HDS, Building, Mesh_facet_handle>;
            using Planar_region_growing = CGAL::LOD::Level_of_detail_planar_region_growing<Kernel, Mesh, Mesh_facets>;

			using All_input_regions = std::vector<Input_regions>;

            using Region_facet   = std::vector<Point_3>;
            using Region_facets  = std::vector<Region_facet>;
            using Output_regions = std::vector<Region_facets>;

            typename Kernel::Compute_squared_length_3         squared_length_3;
            typename Kernel::Compute_scalar_product_3 		  dot_product_3;
			typename Kernel::Construct_cross_product_vector_3 cross_product_3;

            using Color = CGAL::Color;

			Level_of_detail_planar_region_merger() { 

                srand(time(NULL));
            }

			void merge(Clean_facets &clean_facets) const {

                Mesh mesh; Input_regions input_regions;
                create_mesh_and_input_regions(clean_facets, mesh, input_regions);

                Output_regions output_regions;
                create_output_regions(input_regions, output_regions);

                merge_facets(output_regions, clean_facets);
            }

		private:
			
            void create_mesh_and_input_regions(const Clean_facets &clean_facets, Mesh &mesh, Input_regions &input_regions) const {

                mesh.clear();

                Local_builder builder(clean_facets);
                mesh.delegate(builder);

                Mesh_facets mesh_facets;
                builder.get_facets(mesh_facets);

                Planar_region_growing planar_region_growing(mesh_facets);    
                planar_region_growing.use_global_conditions(false);

				All_input_regions all_input_regions;
				planar_region_growing.find_regions(all_input_regions);

                CGAL_postcondition(all_input_regions.size() == 1);
                
                input_regions.clear();
                input_regions = all_input_regions[0];
            }

            void create_output_regions(const Input_regions &input_regions, Output_regions &output_regions) const {
                
                output_regions.clear();
                output_regions.resize(input_regions.size());

                for (size_t i = 0; i < input_regions.size(); ++i)
                    add_output_region(input_regions[i], output_regions[i]);
            }

            void add_output_region(const Input_region &input_region, Region_facets &region_facets) const {
                
                CGAL_precondition(input_region.size() > 0);
                region_facets.clear();

                for (size_t i = 0; i < input_region.size(); ++i)
                    add_output_region_facet(input_region[i], region_facets);
            }

            void add_output_region_facet(const Mesh_facet_handle &fh, Region_facets &region_facets) const {

                Region_facet region_facet;

                Halfedge_handle curr = fh->halfedge();
                Halfedge_handle end  = curr;

                Facet_vertex_handle vh = curr->vertex();
                region_facet.push_back(vh->point());

                curr = curr->next();
                while (curr != end) {

                    vh = curr->vertex();
                    region_facet.push_back(vh->point());

                    curr = curr->next();
                };

                region_facets.push_back(region_facet);
            }

            void merge_facets(Output_regions &output_regions, Clean_facets &clean_facets) const {

                clean_facets.clear();
                for (size_t i = 0; i < output_regions.size(); ++i)
                    merge_region_facets(output_regions[i], clean_facets);
            }

            bool merge_region_facets(Region_facets &region_facets, Clean_facets &clean_facets) const {
                
                const Color color = generate_random_color();

                Vector_3 source_normal;
                compute_source_normal(region_facets, source_normal);

                Vector_3 target_normal;
                compute_target_normal(target_normal);

                if (source_normal == -target_normal) source_normal = target_normal;

                FT angle; Vector_3 axis;
                compute_angle_and_axis(source_normal, target_normal, angle, axis);

                rotate_region_facets(region_facets,  angle, axis);

                CDT cdt;
                triangulate_region_facets(region_facets, cdt);
                get_back_region_facets(cdt, region_facets);
                
                rotate_region_facets(region_facets, -angle, axis);
                
                for (size_t i = 0; i< region_facets.size(); ++i)
                    clean_facets.push_back(std::make_pair(region_facets[i], color));
            }

            void rotate_region_facets(Region_facets &region_facets, const FT angle, const Vector_3 &axis) const {

                for (size_t i = 0; i < region_facets.size(); ++i)
                    rotate_region_facet(region_facets[i], angle, axis);
            }

            void rotate_region_facet(Region_facet &region_facet, const FT angle, const Vector_3 &axis) const {

                Point_3 b;
                compute_barycentre(region_facet, b);

                if (angle != FT(0)) 
                    rotate_facet(b, angle, axis, region_facet);
            }

            void compute_source_normal(const Region_facets &region_facets, Vector_3 &normal) const {
                CGAL_precondition(region_facets.size() > 0);

                Vector_3 tmp_normal; 
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < region_facets.size(); ++i) {
                    compute_facet_normal(region_facets[i], tmp_normal);

                    x += tmp_normal.x();
                    y += tmp_normal.y();
                    z += tmp_normal.z();
                }

                x /= static_cast<FT>(region_facets.size());
                y /= static_cast<FT>(region_facets.size());
                z /= static_cast<FT>(region_facets.size());

                normal = Vector_3(x, y, z);
                normalize(normal);
            }

            void compute_facet_normal(const Region_facet &region_facet, Vector_3 &normal) const {
                CGAL_precondition(region_facet.size() >= 3);
                
                const Point_3 &p1 = region_facet[0];
                const Point_3 &p2 = region_facet[1];
                const Point_3 &p3 = region_facet[2];

                const Vector_3 v1 = Vector_3(p1, p2);
                const Vector_3 v2 = Vector_3(p1, p3);

                normal = cross_product_3(v1, v2);
                normalize(normal);
            }

            void compute_target_normal(Vector_3 &normal) const {
                normal = Vector_3(FT(0), FT(0), FT(1));
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {

				const auto cross = cross_product_3(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length_3(cross))));
				const FT dot     = dot_product_3(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                if (angle == FT(0)) return;

                CGAL_precondition(length != FT(0));
                axis = cross / length;

                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                }
			}

            void compute_barycentre(const Region_facet &region_facet, Point_3 &b) const {

                CGAL_precondition(region_facet.size() != 0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (size_t i = 0; i < region_facet.size(); ++i) {
                    const Point_3 &p = region_facet[i];

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(region_facet.size());
                y /= static_cast<FT>(region_facet.size());
                z /= static_cast<FT>(region_facet.size());

                b = Point_3(x, y, z);
            }

            void rotate_facet(const Point_3 &b, const FT angle, const Vector_3 &axis, Region_facet &region_facet) const {

                Point_3 q;
                for (size_t i = 0; i < region_facet.size(); ++i) {
                    Point_3 &p = region_facet[i];

                    q = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                    rotate_point(angle, axis, q);
                    p = Point_3(q.x() + b.x(), q.y() + b.y(), q.z() + b.z());
                }
            }

            void rotate_point(const FT angle, const Vector_3 &axis, Point_3 &p) const {

				const double tmp_angle = CGAL::to_double(angle);

				const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				const FT C = FT(1) - c;

				const FT x = axis.x();
				const FT y = axis.y();
				const FT z = axis.z();

				p = Point_3((x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
					  		(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
					  		(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
			}

            void triangulate_region_facets(const Region_facets &region_facets, CDT &cdt) const {

            }

            void get_back_region_facets(const CDT &cdt, Region_facets &region_facets) const {

            }

            Color generate_random_color() const {

				const int r = rand() % 255;
				const int g = rand() % 255;
				const int b = rand() % 255;

				return Color(r, g, b);
			}
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PLANAR_REGION_MERGER_H