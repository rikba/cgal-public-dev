#pragma once
#include "defs.h"
#include "defs_cgal.h"

namespace JPTD {

typedef enum {
	INSTANTANEOUS,
	PROGRESSIVE
} Translation_Type;



class Segment_Translation
{
public:
	// Constructor for instantaneous translations
	Segment_Translation(const FT & t, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B);

	// Constructor for progressive translations
	Segment_Translation(const FT & t, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

	~Segment_Translation();

	void set_end(const FT & t);

public:
	const Translation_Type type;
	FT t_int_start;
	FT t_int_end;
	CGAL_Point_2 A;
	CGAL_Point_2 B;
	CGAL_Vector_2 dA;
};

}