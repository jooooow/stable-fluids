#pragma once

#include <vector>

struct Fluids2D
{
	float* u;
	float* v;
	float* p;
	float* u_temp;
	float* v_temp;
	std::vector<float*> scalars;
	std::vector<float*> scalars_temp;
};