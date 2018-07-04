#include "parameters.h"
#include <iostream>
#include <boost/filesystem.hpp>

namespace JPTD {

Kinetic_Parameters::Kinetic_Parameters()
{
	reset();

#ifdef KINETIC_PARTITION_DEBUG
	print_schedule = true;
#endif
}


Kinetic_Parameters::~Kinetic_Parameters()
{
}


void Kinetic_Parameters::reset()
{
	rand_n = 0;
	rand_p = 4;
	rand_d = 0.1;

	path = "";
	location = "";
	basename = "";

	boxes = 10;
	D = FT(1);
	D_inf_2 = FT(9801) / FT(10000); // (0.99 * D)^2
	D_sup_2 = FT(10201) / FT(10000);  // (1.01 * D)^2

	output_facets = false;
	output_polyhedrons = false;

	check = false;
	print_schedule = false;
	print_drawings = false;
}


void Kinetic_Parameters::message() const
{
	std::cout << "Kinetic-Partition-3D.exe [--rand N | --input PATH]" << std::endl;
}



void Kinetic_Parameters::set_options(int argc, char *argv[])
{
	int r = 1;
	while (r < argc) {
		if (!strcmp(argv[r], "--rn") && r + 1 < argc) {
			rand_n = atoi(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--rp") && r + 1 < argc) {
			rand_p = atoi(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--rd") && r + 1 < argc) {
			rand_d = atof(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--input") && r + 1 < argc) {
			path = argv[r + 1];
			location = boost::filesystem::path(path).parent_path().string();
			basename = boost::filesystem::path(path).stem().string();
			r += 2;
		} else if (!strcmp(argv[r], "--box") && r + 1 < argc) {
			boxes = atoi(argv[r + 1]);
			r += 2;
		} else if (!strcmp(argv[r], "--tau") && r + 1 < argc) {
			D = read_quotient(std::string(argv[r + 1]));
			FT D_inf = FT(99) * D / FT(100), D_sup = FT(99) * D / FT(100);
			D_inf_2 = D_inf * D_inf;
			D_sup_2 = D_sup * D_sup;
			r += 2;
		} else if (!strcmp(argv[r], "--polyhedrons")) {
			output_polyhedrons = true;
			r += 1;
		} else if (!strcmp(argv[r], "--facets")) {
			output_facets = true;
			r += 1;
		} else if (!strcmp(argv[r], "--check")) {
			check = true;
			r += 1;
		} else if (!strcmp(argv[r], "--print_schedule")) {
			print_schedule = true;
			r += 1;
		} else if (!strcmp(argv[r], "--print-drawings")) {
			print_drawings = true;
			r += 1;
		} else {
			reset();
			message();
			break;
		}
	}
}



void Kinetic_Parameters::set_option(const std::string & option)
{
	char* args[1] = { strdup(option.c_str()) };
	set_options(1, args);

	free(args[0]);
}



void Kinetic_Parameters::set_option(const std::string & option, const std::string & argument_1)
{
	if (option == "--basename") {
		basename = argument_1;
	} else {
		char* args[2] = { strdup(option.c_str()), strdup(argument_1.c_str()) };
		set_options(2, args);

		free(args[0]);
		free(args[1]);
	}
}



void Kinetic_Parameters::set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2)
{
	char* args[3] = { strdup(option.c_str()), strdup(argument_1.c_str()), strdup(argument_2.c_str()) };
	set_options(3, args);

	free(args[0]);
	free(args[1]);
	free(args[2]);
}



FT Kinetic_Parameters::read_quotient(const std::string & frac)
{	
	FT Q;
	int n = 1, d = 1;
	for (size_t s = 0 ; s < frac.size() ; s++) {
		if (frac[s] == '/') {
			std::string num = frac.substr(0, s);
			std::string den = frac.substr(s + 1, frac.size() - s - 1);
			n = atoi(num.c_str());
			d = atoi(den.c_str());
			Q = FT(n) / FT(d);
			return Q;
		}
	}

	n = atoi(frac.c_str());
	Q = FT(n);
	return Q;
}
}