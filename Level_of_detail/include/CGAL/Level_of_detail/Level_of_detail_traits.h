#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

// STL includes.
#include <map>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// New CGAL includes.
#include <CGAL/Level_of_detail/Level_of_detail_enumerations.h>

#include <CGAL/Level_of_detail/Tools/Tools_include.h>
#include <CGAL/Level_of_detail/Visibility/Visibility_include.h>
#include <CGAL/Level_of_detail/Regularization/Regularization_include.h>
#include <CGAL/Level_of_detail/Shape_detection/Shape_detection_include.h>
#include <CGAL/Level_of_detail/Reconstruction/Reconstruction_include.h>
#include <CGAL/Level_of_detail/Partitioning/Partitioning_include.h>
#include <CGAL/Level_of_detail/Objects/Objects_include.h>

#include "../../../test/Level_of_detail/include/loaders/Myloader.h"
#include "../../../test/Level_of_detail/include/terminal/Myterminal_parser.h"

namespace CGAL {

	namespace Level_of_detail {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits {

			typedef KernelTraits 	Kernel;
			typedef OutputContainer Container_3D;

			typedef CGAL::Level_of_detail::Myloader<Kernel, Container_3D> 					  Loader;
			typedef CGAL::Level_of_detail::Level_of_detail_preprocessor<Kernel, Container_3D> Preprocessor;

			typedef CGAL::Level_of_detail::Level_of_detail_clutter<Kernel, Container_3D> 		   Clutter_strategy;
			typedef CGAL::Level_of_detail::Level_of_detail_ground<Kernel, Container_3D> 		   Ground_strategy;
			typedef CGAL::Level_of_detail::Level_of_detail_building_boundary<Kernel, Container_3D> Building_boundary_strategy;
			typedef CGAL::Level_of_detail::Level_of_detail_building_interior<Kernel, Container_3D> Building_interior_strategy;

			typedef CGAL::Level_of_detail::Level_of_detail_selector<Kernel, Clutter_strategy> 		    Clutter_selector;
			typedef CGAL::Level_of_detail::Level_of_detail_selector<Kernel, Ground_strategy> 		    Ground_selector;
			typedef CGAL::Level_of_detail::Level_of_detail_selector<Kernel, Building_boundary_strategy> Building_boundary_selector;
			typedef CGAL::Level_of_detail::Level_of_detail_selector<Kernel, Building_interior_strategy> Building_interior_selector;

			typedef std::map<int, std::vector<int> >          				   		   			  		Planes;
			typedef std::map<int, typename Kernel::Point_2>           	    				  		 	Projected_points;
			typedef CGAL::Level_of_detail::Level_of_detail_simple_projector<Kernel, Container_3D, Projected_points> Ground_projector;
			
			typedef CGAL::Level_of_detail::Level_of_detail_segment_regularizer_2<Kernel> Line_regularizer;

			typedef CGAL::Level_of_detail::My_vertex_info<Structured_label>  My_vertex_info; 
	    	typedef CGAL::Level_of_detail::My_face_info<typename Kernel::FT> My_face_info;

			typedef CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Kernel> VB;
	 	    typedef CGAL::Triangulation_face_base_with_info_2<My_face_info, Kernel>     FB_with_info;
			typedef CGAL::Constrained_triangulation_face_base_2<Kernel, FB_with_info>   FB;

			   typedef CGAL::Exact_predicates_tag                                TAG;
			// typedef CGAL::Exact_intersections_tag                             TAG;

			typedef CGAL::Triangulation_data_structure_2<VB, FB> 			  	 TDS;
			typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, TAG> CDT;
			
			// typedef CGAL::Constrained_triangulation_plus_2<DDT> CDT;

			typedef CGAL::Level_of_detail::Level_of_detail_utils<Kernel, Container_3D, CDT> 				   				 Utils;
			typedef CGAL::Level_of_detail::Level_of_detail_region_growing_2<Kernel, Planes, Projected_points, Container_3D>  Region_growing_2;

			typedef CGAL::Level_of_detail::Level_of_detail_grid_simplify<Kernel, Planes, Projected_points> Grid_simplifier;

			typedef int Label;

			typedef std::pair<typename Kernel::Point_2, Label> Point_with_label;
			typedef std::vector<Point_with_label> 			   Container_2D;
			
			typedef CGAL::Level_of_detail::Building<Kernel, typename CDT::Vertex_handle, typename CDT::Face_handle> Building;
			typedef std::map<int, Building> 															            Buildings;

			typedef CGAL::Polyhedron_3<Kernel> 											                Mesh;
			typedef CGAL::Level_of_detail::Level_of_detail_reconstruction<Kernel, CDT, Buildings, Mesh> Lods;

			typedef typename Lods::Mesh_facet_colors Mesh_facet_colors;

			typedef CGAL::Level_of_detail::Level_of_detail_building_splitter<Kernel, CDT> Building_splitter;
			typedef CGAL::Level_of_detail::Level_of_detail_building_outliner<Kernel, CDT> Building_outliner;

			typedef CGAL::Level_of_detail::Level_of_detail_min_height_fitter<Kernel> Min_height_fitter;
			typedef CGAL::Level_of_detail::Level_of_detail_avg_height_fitter<Kernel> Avg_height_fitter;
			typedef CGAL::Level_of_detail::Level_of_detail_max_height_fitter<Kernel> Max_height_fitter;
	
			typedef CGAL::Level_of_detail::Level_of_detail_building_roof_fitter<Kernel, CDT, Container_3D, Min_height_fitter> Building_min_roof_fitter;
			typedef CGAL::Level_of_detail::Level_of_detail_building_roof_fitter<Kernel, CDT, Container_3D, Avg_height_fitter> Building_avg_roof_fitter;
			typedef CGAL::Level_of_detail::Level_of_detail_building_roof_fitter<Kernel, CDT, Container_3D, Max_height_fitter> Building_max_roof_fitter;

			typedef CGAL::Level_of_detail::Myterminal_parser<typename Kernel::FT> Level_of_detail_parameters;
			typedef typename Level_of_detail_parameters::Input_parameters 	      Parameters;

			typedef CGAL::Level_of_detail::Level_of_detail_parameters_estimator<Kernel, Container_3D, Parameters> Parameters_estimator;

			typedef CGAL::Level_of_detail::Level_of_detail_complexity<Kernel, Container_3D, Lods> Lod_complexity;
			typedef CGAL::Level_of_detail::Level_of_detail_distortion<Kernel, Container_3D, Lods> Lod_distortion;
			typedef CGAL::Level_of_detail::Level_of_detail_coverage<Kernel, Container_3D, Lods>   Lod_coverage;

			typedef CGAL::Level_of_detail::Level_of_detail_container<Kernel> 								     Lod_data_structure;
			typedef CGAL::Level_of_detail::Level_of_detail_polygonizer_jean_philippe<Kernel, Lod_data_structure> Polygonizer;
			
            typedef std::pair<typename Kernel::FT, typename Kernel::FT> Visibility_pair;
            typedef std::map<size_t, Visibility_pair> 					Visibility_output;

			typedef CGAL::Level_of_detail::Level_of_detail_classification_shepard_visibility_strategy_2<Kernel, Container_2D, Lod_data_structure, Visibility_output> Visibility_strategy;
			typedef CGAL::Level_of_detail::Level_of_detail_polygon_based_visibility_2<Kernel, Container_3D, Lod_data_structure, Visibility_strategy> 			     Polygon_based_visibility;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H