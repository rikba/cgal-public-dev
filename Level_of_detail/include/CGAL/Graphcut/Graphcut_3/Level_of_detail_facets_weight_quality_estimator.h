#ifndef CGAL_LEVEL_OF_DETAIL_FACETS_WEIGHT_QUALITY_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_FACETS_WEIGHT_QUALITY_ESTIMATOR_H 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class InputKernel, class InputData, class InputBuilding>
		class Level_of_detail_facets_weight_quality_estimator {

		public:
            using Kernel   = InputKernel;
			using Input    = InputData;
			using Building = InputBuilding;

			using FT 	  = typename Kernel::FT;
			using Point_3 = typename Kernel::Point_3;

			Level_of_detail_facets_weight_quality_estimator(const Input &input, const FT ground_height) : 
			m_input(input),
			m_ground_height(ground_height) { }

			void estimate(Building &building) const {

			}

		private:
			const Input &m_input;
			const FT m_ground_height;
		};

	} // CGAL

} // LOD

#endif // CGAL_LEVEL_OF_DETAIL_FACETS_WEIGHT_QUALITY_ESTIMATOR_H