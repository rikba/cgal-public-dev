#ifndef CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H
#define CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H

#if defined(WIN32) || defined(_WIN32) 
#define _SR_ "\\"
#else 
#define _SR_ "/" 
#endif

// STL includes.
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail_include.h>

// Local includes.
#include "debugging/Mylog.h"
#include "loaders/Myloader.h"
#include "terminal/Myterminal_parser.h"

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel>
		class Mywrapper {

		public:
			using Kernel = InputKernel;
			
			using Parameters = char**;
			using FT 		 = typename Kernel::FT;
			
			using Log 			  = LOD::Mylog;
			using Loader 		  = LOD::Myloader<Kernel>;
			using Terminal_parser = LOD::Myterminal_parser<FT>;

			using Label = typename Loader::Label;
			
			using LOD_input     = typename Loader::Container;
			using LOD_point_map = typename Loader::Point_map;
			using LOD_label_map = typename Loader::Label_map;

			using LOD_point 	  = typename Loader::Point_3;
			using LOD_point_index = typename Loader::Point_index;

			using LOD_semantic_element_map = LOD::Semantic_element_property_map<LOD_point_index, LOD_label_map>;

			using LOD_base 		 = LOD::Level_of_detail<Kernel, LOD_input, LOD_point_map>;
			using LOD_parameters = typename LOD_base::Parameters;

			using LOD_dereference_point_map = LOD::Dereference_property_map<typename LOD_base::Point_identifier, LOD_point_map>;

			Mywrapper(const int num_parameters, const Parameters parameters, const std::string &logs_path) : 
			m_terminal_parser(num_parameters, parameters, logs_path),
			m_logs_path_0_1(logs_path + "tmp" + std::string(_SR_) + "lod_0_1" + std::string(_SR_))
			{ }

			void run_lod_pipeline() {
				
				// First, we set LOD parameters parsing a terminal input or to the default values, including the path to the input data.
				set_lod_parameters();

				// Then, we load data from a file. I use my own custom loader here.
				load_lod_input_data();

				// We also set all necessary property maps.
				set_lod_property_maps();

				// Here, we build all steps of the LOD pipeline with the default template classes.
				// compute_lod_default();
				
				// Here, we build all steps of the LOD pipeline with the user-defined template classes.
				compute_lod_custom();
			}

		private:
			Terminal_parser m_terminal_parser;
			LOD_parameters  m_lod_parameters;
			LOD_input       m_lod_input;

			LOD_point_map m_lod_point_map;
			LOD_label_map m_lod_label_map;

			std::string m_logs_path_0_1;

			void set_lod_parameters() {

				// Set all parameters that can be loaded from the terminal.
				// add_str_parameter  - adds a string-type parameter
				// add_val_parameter  - adds a scalar-type parameter
				// add_bool_parameter - adds a boolean parameter

				std::cout << "Input parameters: " << std::endl;

				m_terminal_parser.add_str_parameter("-data", m_lod_parameters.path_to_input());
				m_terminal_parser.add_bool_parameter("-verbose", m_lod_parameters.verbose());
				
				m_terminal_parser.add_val_parameter("-scale", m_lod_parameters.scale());
				m_lod_parameters.update_scale_dependent();
			}

			void load_lod_input_data() {

				std::cout << std::endl << "Input data: " << std::endl;
				Loader loader; loader.get_data(m_lod_parameters.path_to_input(), m_lod_input);

				std::cout << "LOD input data are loaded: number of points: " << m_lod_input.number_of_points() << std::endl;
			}

			void set_lod_property_maps() {

				m_lod_point_map = m_lod_input.point_map();
				m_lod_label_map = m_lod_input. template property_map<Label>("label").first;
			}

			void compute_lod_default() {
				
				// Begin.
				CGAL::Timer timer;

				std::cout << std::endl << "Default LOD pipeline: " << std::endl;
				LOD_base lod_base(m_lod_input, m_lod_point_map, m_lod_parameters);
				
				timer.start();

				LOD_semantic_element_map lod_semantic_element_map(m_lod_label_map);
				lod_base.build(lod_semantic_element_map);

				std::cout << std::endl;
				lod_base.get_lod0();
				lod_base.get_lod1();

				// End.
				timer.stop();
				std::cout << std::endl << "Running time: " << timer.time() << " seconds." << std::endl << std::endl;
			}

			void compute_lod_custom() {

				// Begin.
				Log log;
				CGAL::Timer timer;

				std::cout << std::endl << "Custom LOD pipeline: " << std::endl;
				LOD_base lod_base(m_lod_input, m_lod_point_map, m_lod_parameters);

				timer.start();

				// * Step ->
				LOD_semantic_element_map lod_semantic_element_map(m_lod_label_map);
				lod_base.split_semantic_data(lod_semantic_element_map);
				
				LOD_dereference_point_map lod_dereference_point_map(m_lod_point_map);

				log.save_points(lod_base.get_internal_data_structure().ground_points(), lod_dereference_point_map, m_logs_path_0_1 + "0_ground_points");
				log.save_points(lod_base.get_internal_data_structure().building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "1_building_boundary_points");
				log.save_points(lod_base.get_internal_data_structure().building_interior_points(), lod_dereference_point_map, m_logs_path_0_1 + "2_building_interior_points");

				// * Step ->
				lod_base.fit_ground_plane();

				Identity_property_map<LOD_point> point_map;
				log.save_polygon(lod_base.get_internal_data_structure().ground_bounding_box(), point_map, m_logs_path_0_1 + "3_ground_plane");

				// * Step ->
				lod_base.extract_building_boundaries();
				log.save_points(lod_base.get_internal_data_structure().filtered_building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "4_filtered_building_boundary_points");

				// * Step ->
				lod_base.simplify_building_boundaries();
				log.save_points(lod_base.get_internal_data_structure().simplified_building_boundary_points(), lod_dereference_point_map, m_logs_path_0_1 + "5_simplified_building_boundary_points");

				// * Step ->
				std::cout << std::endl;
				lod_base.get_lod0();
				lod_base.get_lod1();

				// End.
				timer.stop();
				std::cout << std::endl << "Running time: " << timer.time() << " seconds." << std::endl << std::endl;
			}
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MYWRAPPER_H