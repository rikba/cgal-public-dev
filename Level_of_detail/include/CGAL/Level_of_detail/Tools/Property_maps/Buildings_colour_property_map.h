#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_COLOUR_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_COLOUR_PROPERTY_MAP_H

// STL includes.
#include <map>

// CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Property_maps/Colour_property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		class Buildings_colour_property_map {

		public:
            using Colour    = CGAL::Color;
            using ValueType = Colour;

            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Buildings_colour_property_map;

            using Building_colours          = std::map<int, Colour>;
            using Building_colours_iterator = typename Building_colours::const_iterator;

            using Colour_map      = LOD::Colour_property_map;
            using Colour_map_type = LOD::Colour_map_type;

            Buildings_colour_property_map(const size_t num_buildings) {

                for (size_t i = 0; i < num_buildings; ++i)
                    m_building_colours[i] = get_new_colour(i);
            }

            inline const Building_colours &building_colours() const {
                return m_building_colours;
            }

			template<typename key_type>
			value_type operator[](key_type &key) const { 
				return get(this, key);
			}

			template<typename key_type>
            friend value_type get(const Self &self, const key_type &key) {
				return self.generate_building_colour(self, key);
            }

			template<class Facet>
			Colour generate_building_colour(const Self &self, const Facet &facet) const {
				
                const int building_number = facet.info().group_number();
                if (building_number < 0) return self.get_default_colour(building_number);

                CGAL_precondition(does_building_colour_exist(building_number, self.building_colours()));
                return self.building_colours().at(building_number);
			}

        private:
            Building_colours m_building_colours;

            Colour get_default_colour(const int building_number) const {
                
                const Colour_property_map colour_map(Colour_map_type::WHITE);
                return get(colour_map, building_number);
            }

            Colour get_new_colour(const int building_number) const {
                
                const Colour_property_map colour_map(Colour_map_type::RANDOM);
                return get(colour_map, building_number);
            }

            inline bool does_building_colour_exist(const int building_number, const Building_colours &building_colours) const {
                return building_number < building_colours.size();
            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_COLOUR_PROPERTY_MAP_H