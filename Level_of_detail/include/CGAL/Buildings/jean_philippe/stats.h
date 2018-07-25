#pragma once
#include <ctime>

namespace JPTD {

namespace KP_Stats 
{
	extern int insert_temporally_calls;
	extern int insert_temporally_loop_case;
	extern int insert_temporally_loop_comparisons;

	extern clock_t schedule_events_search_time;
	extern clock_t schedule_events_computation_time;

	extern int schedule_events_vertices;
	extern int schedule_events_lines;

	extern int life_expectancy_lines;
	extern double life_expectancy_distance;
	extern int life_expectancy_terms;

	extern clock_t process_events_time;
	extern int process_events_calls;
}

}