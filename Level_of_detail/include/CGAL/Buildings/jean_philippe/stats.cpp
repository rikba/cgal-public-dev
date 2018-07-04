#include "stats.h"

namespace JPTD {
int KP_Stats::insert_temporally_calls = 0;
int KP_Stats::insert_temporally_loop_case = 0;
int KP_Stats::insert_temporally_loop_comparisons = 0;

clock_t KP_Stats::schedule_events_search_time = 0;
clock_t KP_Stats::schedule_events_computation_time = 0;

int KP_Stats::schedule_events_vertices = 0;
int KP_Stats::schedule_events_lines = 0;

int KP_Stats::life_expectancy_lines = 0;
double KP_Stats::life_expectancy_distance = 0;
int KP_Stats::life_expectancy_terms = 0;

clock_t KP_Stats::process_events_time = 0;
int KP_Stats::process_events_calls = 0;
}