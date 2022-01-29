#pragma once

#include "type.h"

void UpdateFluids(
	int width,
	int height,
	unsigned char* region_flag,
	Fluids2D* fluids
);

void AdvectVelocity(
	int width,
	int height,
	unsigned char* region_flag,
	float* u,
	float* v,
	float* u_temp,
	float* v_temp
);

void ProjectVelocity(
	int width,
	int height,
	unsigned char* region_flag,
	float* u,
	float* v,
	float* u_temp,
	float* v_temp,
	float* p
);

void AdvectScalar(
	int width,
	int height,
	unsigned char* region_flag,
	float* u,
	float* v,
	float* scalar,
	float* scalar_temp
);

void BoundaryConditionVelocity(
	int width,
	int height,
	unsigned char* region_flag,
	float* u,
	float* v
);

float Inter2(float t, float s, float a, float b, float c, float d);