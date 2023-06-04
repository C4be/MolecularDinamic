#pragma once
#include "global_vars.h"

void initial_condition_two_particles()
{
	coordx[0] = 0.75; coordy[0] = 0.75; coordz[0] = 0.5;
	coordx[1] = 1.25; coordy[1] = 0.75; coordz[1] = 0.5;

	vx[0] = 1.0; vy[0] = 1.0; vz[0] = 0.0;
	vx[1] = -1.0; vy[1] = 1.0; vz[1] = 0.0;
}