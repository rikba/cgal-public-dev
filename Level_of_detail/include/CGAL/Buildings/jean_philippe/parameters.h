#pragma once
#include "defs.h"
#include <string>

namespace JPTD {
class Kinetic_Parameters
{
public:
	Kinetic_Parameters();

	~Kinetic_Parameters();

	void message() const;

	void set_options(int argc, char *argv[]);

	void set_option(const std::string & option);

	void set_option(const std::string & option, const std::string & argument_1);

	void set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2);

protected:
	void reset();

	FT read_quotient(const std::string & frac);

public:
	int rand_n;
	int rand_p;
	double rand_d;

	std::string path;
	std::string location;
	std::string basename;

	uint boxes;
	FT D;
	FT D_inf_2;
	FT D_sup_2;

	bool output_facets;
	bool output_polyhedrons;

	bool check;
	bool print_schedule;
	bool print_drawings;
};
}