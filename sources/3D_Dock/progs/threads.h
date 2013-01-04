
#include <pthread.h>
#include <semaphore.h>

#ifndef __DEF_THREADS_
#define __DEF_THREADS_

#ifndef NTHREADS
#define NTHREADS 6
#endif

extern sem_t num_threads_sem;

struct rfftw3d_create_plan_parameters {
	rfftwnd_plan *ret;
	int nx;
	int ny;
	int nz;
	fftw_direction dir;
	int flags;
};

struct electric_field_parameters {
	struct Structure This_Structure;
	float grid_span;
	int grid_size;
	fftw_real *grid;
	int *shared_x;
	struct atom_values *atoms;
	int natoms_in;
};

void *rfftw3d_create_plan_thread(void *v_params);
void *electric_field_thread(void *v_params);

#endif

