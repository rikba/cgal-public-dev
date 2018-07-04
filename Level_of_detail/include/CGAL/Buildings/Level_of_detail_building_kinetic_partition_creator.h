#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

// Jean Philippe includes.
#include <CGAL/Buildings/jean_philippe/propagation.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuildings, class InputBuilding>
		class Level_of_detail_building_kinetic_partition_creator {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Segment_2  = typename Kernel::Segment_2;
            using Segment_3  = typename Kernel::Segment_3;

            using CDT   = typename Building::CDT;
            using Roofs = typename Building::Roofs;
            
            using Roof = typename Building::Roof;
            using Roof_boundary = typename Roof::Roof_boundary;

            using Building_iterator = typename Buildings::iterator;
            using Building_boundary = typename Building::Boundary;
            
            using Log      = CGAL::LOD::Mylog;
            using Segments = std::vector<Segment_2>;

            using Polygon_boundary = std::vector<Point_3>;
            using Polygons = std::vector<Polygon_boundary>;

            Level_of_detail_building_kinetic_partition_creator(const CDT &cdt, const FT ground_height) :
            m_cdt(cdt),
            m_ground_height(ground_height)
            { }

            void create(Buildings &buildings) const {
                
                if (buildings.size() == 0) return; int count = 0;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit, ++count) {

                    Building &building = bit->second; building.index = count;
					if (building.is_valid) process_building(building);
                }
            }

        private:
            const CDT &m_cdt;
            const FT m_ground_height;
            
            void process_building(Building &building) const {
                
                const std::string path = save_convex_polygons(building);
                std::cout << std::endl << "input file: " << path << std::endl << std::endl;

                const std::string str1 = "./kinetic";
                const std::string str2 = "--input";
                const std::string str3 = "--facets";

                char *val1 = strdup(str1.c_str());
                char *val2 = strdup(str2.c_str());
                char *val3 = strdup(path.c_str());
                char *val4 = strdup(str3.c_str());

                int argc = 4;
                char *n_argv[] = { val1, val2, val3, val4 };
                char **argv = n_argv;

                JPTD::Kinetic_Propagation kinetic(argc, argv);
	            if (!kinetic.data()) return;
                kinetic.run();
            }

            std::string save_convex_polygons(const Building &building) const {

                Polygons polygons;
                const Building_boundary &building_boundary = building.boundaries[0];

                // Add boundaries.
                for (size_t i = 0; i < building_boundary.size(); i += 2) {
					const size_t ip = i + 1;

					assert(ip < building_boundary.size());
                    const Point_2 &a = building_boundary[i]->point();
                    const Point_2 &b = building_boundary[ip]->point();

                    const Point_3 p1 = Point_3(a.x(), a.y(), m_ground_height);
                    const Point_3 p2 = Point_3(b.x(), b.y(), m_ground_height);
                    const Point_3 p3 = Point_3(b.x(), b.y(), m_ground_height + FT(2));
                    const Point_3 p4 = Point_3(a.x(), a.y(), m_ground_height + FT(2));

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

                Log log;

                const std::string path = "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "kinetic_input_building_" + std::to_string(building.index);
                log.save_convex_polygons(polygons, path);

                return log.m_prefix_path + path + ".ply";
            }
        };

    } // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_KINETIC_PARTITION_CREATOR_H