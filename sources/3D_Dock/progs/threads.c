
#include "structures.h"
#include "threads.h"

void *rfftw3d_create_plan_thread(void *v_params)
{
	struct rfftw3d_create_plan_parameters *params = v_params;

	*params->ret = rfftw3d_create_plan(params->nx, params->ny, params->nz, params->dir, params->flags);
	
	return NULL;
}

