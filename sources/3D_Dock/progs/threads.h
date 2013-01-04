
#include <pthread.h>

#ifndef __DEF_THREADS_
#define __DEF_THREADS_

#ifndef NTHREADS
#define NTHREADS 4
#endif

void *rfftw3d_create_plan_thread(void *v_params);

#endif

