#include "fluids.h"
#include <math.h>
#include <iostream>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define dt 2
#define dx 1
#define rou 1
#define RELAX_ITER 300

void UpdateFluids(int width, int height, unsigned char* region_flag, Fluids2D* fluids)
{
	// update velocity field
	AdvectVelocity(width, height, region_flag, fluids->u, fluids->v, fluids->u_temp, fluids->v_temp);  // new u, v are stored in u_temp, v_temp
	BoundaryConditionVelocity(width, height, region_flag, fluids->u_temp, fluids->v_temp);
	memcpy(fluids->u, fluids->u_temp, (width + 1) * height * sizeof(float));
	memcpy(fluids->v, fluids->v_temp, width * (height + 1) * sizeof(float));
	ProjectVelocity(width, height, region_flag, fluids->u, fluids->v, fluids->u_temp, fluids->v_temp, fluids->p);

	// update scalar field
	for (int i = 0; i < fluids->scalars.size(); i++)
	{
		AdvectScalar(width, height, region_flag, fluids->u, fluids->v, fluids->scalars.at(i), fluids->scalars_temp.at(i));
		memcpy(fluids->scalars.at(i), fluids->scalars_temp.at(i), width * height * sizeof(float));
	}
}

void AdvectVelocity(int width, int height, unsigned char* region_flag, float* u, float* v, float* u_temp, float* v_temp)
{
	//std::cout << "AdvectVelocity" << std::endl;
	int uwidth = width + 1;
	int uheight = height;
	int vwidth = width;
	int vheight = height + 1;

	// advect u
	for (int j = 0; j < uheight; j++)
	{
		for (int i = 0; i < uwidth; i++)
		{
			float u_here, v_here;

			u_here = u[j * uwidth + i];
			if (i == 0 || i == uwidth - 1)
			{
				v_here = 0.0f;
			}
			else
			{
				v_here = 0.25f * (v[j * vwidth + i - 1] + v[j * vwidth + i] + v[(j + 1) * vwidth + i - 1] + v[(j + 1) * vwidth + i]);
			}

			float x = i * dx;
			float y = j * dx + 0.5f * dx;
			float x_prev = x - dt * u_here;
			float y_prev = y - dt * v_here;

			if (x_prev < 0 || x_prev >= width * dx || y_prev < 0.5f * dx || y_prev >= height * dx - 0.5 * dx)  // clamp y half cell away from the border since the u-node is in the center of a vertical cell border
			{
				u_temp[j * uwidth + i] = 0.0f;
			}
			else
			{
				int i0 = (int)(x_prev / dx);
				int i1 = i0 + 1;
				int j0 = (int)((y_prev - 0.5f * dx) / dx);
				int j1 = j0 + 1;

				float t = (x_prev - i0 * dx) / dx;
				float s = (y_prev - (j0 * dx + 0.5f * dx)) / dx;

				float a = u[j0 * uwidth + i0];
				float b = u[j0 * uwidth + i1];
				float c = u[j1 * uwidth + i0];
				float d = u[j1 * uwidth + i1];
				u_temp[j * uwidth + i] = Inter2(t, s, a, b, c, d);
			}
		}
	}

	// advect v
	for (int j = 0; j < vheight; j++)
	{
		for (int i = 0; i < vwidth; i++)
		{
			float v_here, u_here;

			v_here = v[j * vwidth + i];

			if (j == 0 || j == vheight - 1)
			{
				u_here = 0.0f;
			}
			else
			{
				u_here = 0.25f * (u[(j - 1) * uwidth + i] + u[(j - 1) * uwidth + i + 1] + u[j * uwidth + i] + u[j * uwidth + i + 1]);
			}

			float x = i * dx + 0.5f * dx;
			float y = j * dx;
			float x_prev = x - dt * u_here;
			float y_prev = y - dt * v_here;

			if (x_prev < 0.5f * dx || x_prev >= width * dx - 0.5f * dx || y_prev < 0 || y_prev >= height * dx)  // clamp x half cell away from the border since the v-node is in the center of a horizontal cell border
			{
				v_temp[j * vwidth + i] = 0.0f;
			}
			else
			{
				int i0 = (int)((x_prev - 0.5f * dx) / dx);
				int i1 = i0 + 1;
				int j0 = (int)(y_prev / dx);
				int j1 = j0 + 1;

				float t = (x_prev - (i0 * dx + 0.5f * dx)) / dx;
				float s = (y_prev - j0 * dx) / dx;

				float a = v[j0 * vwidth + i0];
				float b = v[j0 * vwidth + i1];
				float c = v[j1 * vwidth + i0];
				float d = v[j1 * vwidth + i1];
				v_temp[j * vwidth + i] = Inter2(t, s, a, b, c, d);
			}
		}
	}
}

void ProjectVelocity(int width, int height, unsigned char* region_flag, float* u, float* v, float* u_temp, float* v_temp, float* p)
{
	//std::cout << "ProjectVelocity" << std::endl;
	int uwidth = width + 1;
	int uheight = height;
	int vwidth = width;
	int vheight = height + 1;

	memset(p, 0, width * height * sizeof(float));

	

	// solve pressure
	for (int iter = 0; iter < RELAX_ITER; iter++)
	{
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				if (region_flag[j * width + i] != 0)
				{
					float u1, u2, v1, v2;
					int pij = -4;
					float p_sum = 0.0f;

					if (region_flag[(j - 1) * width + i] == 0)  //p(i, j - 1)
					{
						v1 = 0;  //v_solid
						pij++;
					}
					else
					{
						v1 = v[j * vwidth + i];
						p_sum += p[(j - 1) * width + i];
					}

					if (region_flag[j * width + i + 1] == 0)    //p(i + 1, j)
					{
						u2 = 0;  //u_solid
						pij++;
					}
					else
					{
						u2 = u[j * uwidth + i + 1];
						p_sum += p[j * width + i + 1];
					}

					if (region_flag[(j + 1) * width + i] == 0)  //p(i, j + 1)
					{
						v2 = 0;  //v_solid
						pij++;
					}
					{
						v2 = v[(j + 1) * vwidth + i];
						p_sum += p[(j + 1) * width + i];
					}

					if (region_flag[j * width + i - 1] == 0)    //p(i - 1, j)
					{
						u1 = 0;  //u_soild
						pij++;
					}
					else
					{
						u1 = u[j * uwidth + i];
						p_sum += p[j * width + i - 1];
					}

					float divergence = (u2 - u1 + v2 - v1) * rou * dx / dt;

					if (pij != 0)
					{
						p[j * width + i] = (divergence - p_sum) / pij;
					}
					else
					{
						p[j * width + i] = 0;
					}
				}
			}
		}	
	}

	float resdiual = 0;
	int cnt = 0;

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if (region_flag[j * width + i] != 0)
			{
				cnt++;

				float u1, u2, v1, v2;
				int pij = -4;
				float p_sum = 0.0f;

				if (region_flag[(j - 1) * width + i] == 0)  //p(i, j - 1)
				{
					v1 = 0;  //v_solid
					pij++;
				}
				else
				{
					v1 = v[j * vwidth + i];
					p_sum += p[(j - 1) * width + i];
				}

				if (region_flag[j * width + i + 1] == 0)    //p(i + 1, j)
				{
					u2 = 0;  //u_solid
					pij++;
				}
				else
				{
					u2 = u[j * uwidth + i + 1];
					p_sum += p[j * width + i + 1];
				}

				if (region_flag[(j + 1) * width + i] == 0)  //p(i, j + 1)
				{
					v2 = 0;  //v_solid
					pij++;
				}
				{
					v2 = v[(j + 1) * vwidth + i];
					p_sum += p[(j + 1) * width + i];
				}

				if (region_flag[j * width + i - 1] == 0)    //p(i - 1, j)
				{
					u1 = 0;  //u_soild
					pij++;
				}
				else
				{
					u1 = u[j * uwidth + i];
					p_sum += p[j * width + i - 1];
				}

				float divergence = (u2 - u1 + v2 - v1) * rou * dx / dt;

				resdiual = MAX(divergence - (p_sum + p[j * width + i] * pij), resdiual);
			}
		}
	}

	//resdiual /= cnt;

	printf("resdiual : % .10f\n", resdiual);
	
	

	// subtract divergence
	for (int j = 0; j < uheight; j++)
	{
		for (int i = 0; i < uwidth; i++)
		{
			if (i == 0 || i == uwidth - 1)
			{
				u[j * uwidth + i] = 0;
			}
			else
			{
				if (region_flag[j * width + i - 1] == 0 || region_flag[j * width + i] == 0)
				{
					u[j * uwidth + i] = 0;
				}
				else
				{
					u[j * uwidth + i] -= dt * (p[j * width + i] - p[j * width + i - 1]) / dx / rou;
				}
			}
		}
	}

	for (int j = 0; j < vheight; j++)
	{
		for (int i = 0; i < vwidth; i++)
		{
			if (j == 0 || j == vheight - 1)
			{
				v[j * vwidth + i] = 0;
			}
			else
			{
				if (region_flag[(j - 1) * vwidth + i] == 0 || region_flag[j * vwidth + i] == 0)
				{
					v[j * vwidth + i] = 0;
				}
				else
				{
					v[j * vwidth + i] -= dt * (p[j * width + i] - p[(j - 1) * width + i]) / dx / rou;
				}
			}
		}
	}
}

void AdvectScalar(int width, int height, unsigned char* region_flag, float* u, float* v, float* scalar, float* scalar_temp)
{
	//std::cout << "AdvectScalar" << std::endl;
	int uwidth = width + 1;
	int uheight = height;
	int vwidth = width;
	int vheight = height + 1;

	// advect
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if (region_flag[j * width + i] != 0)
			{
				float u_here, v_here;

				u_here = 0.5f * (u[j * uwidth + i] + u[j * uwidth + i + 1]);
				v_here = 0.5f * (v[j * vwidth + i] + v[(j + 1) * vwidth + i]);

				float x = i * dx + 0.5f * dx;
				float y = j * dx + 0.5f * dx;
				float x_prev = x - dt * u_here;
				float y_prev = y - dt * v_here;

				if (x_prev <0.5f * dx || x_prev > width * dx - 0.5f * dx || y_prev < 0.5f * dx || y_prev > height * dx - 0.5f * dx)
				{
					scalar_temp[j * width + i] = 0.0f;
				}
				else
				{
					int i0 = (int)((x_prev - 0.5f * dx) / dx);
					int i1 = i0 + 1;
					int j0 = (int)((y_prev - 0.5f * dx) / dx);
					int j1 = j0 + 1;

					float t = (x_prev - (i0 * dx + 0.5f * dx)) / dx;
					float s = (y_prev - (j0 * dx + 0.5f * dx)) / dx;

					scalar_temp[j * width + i] = Inter2(t, s, scalar[j0 * width + i0], scalar[j0 * width + i1], scalar[j1 * width + i0], scalar[j1 * width + i1]);
				}
			}
		}
	}
}

void BoundaryConditionVelocity(
	int width,
	int height,
	unsigned char* region_flag,
	float* u,
	float* v
)
{
	int uwidth = width + 1;
	int uheight = height;
	int vwidth = width;
	int vheight = height + 1;

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			if (region_flag[j * width + i] == 0)
			{
				u[j * uwidth + i] = 0.0f;
				u[j * uwidth + i + 1] = 0.0f;
				v[j * vwidth + i] = 0.0f;
				v[(j + 1) * vwidth + i] = 0.0f;
			}
		}
	}
}

float Inter2(float t, float s, float a, float b, float c, float d)
{
	float x0 = (1 - t) * a + t * b;
	float x1 = (1 - t) * c + t * d;
	float y = (1 - s) * x0 + s * x1;
	return y;
}