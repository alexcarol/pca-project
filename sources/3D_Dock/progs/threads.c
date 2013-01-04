
#include "structures.h"
#include "threads.h"

sem_t num_threads_sem;

void *rfftw3d_create_plan_thread(void *v_params)
{
	struct rfftw3d_create_plan_parameters *params = v_params;

	*params->ret = rfftw3d_create_plan(params->nx, params->ny, params->nz, params->dir, params->flags);
	
	sem_post(&num_threads_sem);
	
	return NULL;
}

void *electric_field_thread(void *v_params)
{
	struct electric_field_parameters *params = v_params;

	electric_field(params->This_Structure, params->grid_span, params->grid_size, params->grid, params->shared_x, params->atoms, params->natoms_in);
	
	sem_post(&num_threads_sem);
	
	return NULL;
}

