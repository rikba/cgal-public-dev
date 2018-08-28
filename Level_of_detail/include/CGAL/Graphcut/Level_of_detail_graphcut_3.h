#ifndef CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H
#define CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/IO/Color.h>

namespace Maxflow {
	#include <CGAL/internal/auxiliary/graph.h>
}

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding, class InputBuildings>
		class Level_of_detail_graphcut_3 {

		public:
            using Kernel    = InputKernel;
			using Input     = InputData;
			using Building  = InputBuilding;
			using Buildings = InputBuildings;

			typedef typename Kernel::FT 	 FT;
			typedef typename Kernel::Point_3 Point_3;

			typedef Maxflow::Graph 			Graph;
			typedef typename Graph::node_id Node_id;

			using Buildings_iterator = typename Buildings::iterator;

			using Log = CGAL::LOD::Mylog;

			Level_of_detail_graphcut_3(const Input &input, const FT ground_height) : 
			m_input(input),
			m_ground_height(ground_height),
			m_alpha(FT(1)),
			m_beta(FT(100000)),
			m_gamma(FT(1000)) { 

				
			}

			void set_alpha(const FT alpha) {
				// m_alpha = alpha;
			}

			void set_beta(const FT beta) {
				// m_beta = beta;
			}

			void set_gamma(const FT gamma) {
				// m_gamma = gamma;
			}

			void solve(Buildings &buildings) const {

				if (buildings.size() == 0)
                    return;
                    
				for (Buildings_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    Building &building = bit->second;

					if (building.is_valid && building.shapes.size() != 0 && building.polyhedrons.size() != 0) 
						process_building(building);
					else 
						building.is_valid = false;
                }
			}

		private:
			const Input &m_input;
			const FT m_ground_height;

			FT m_alpha;
			FT m_beta;
			FT m_gamma;

			void process_building(Building &building) const {

                
            }
		};

	} // CGAL

} // LOD

#endif // CGAL_LEVEL_OF_DETAIL_GRAPHCUT_3_H