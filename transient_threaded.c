/* initial conditions (shape is square):
 *
 *       |<-- 1m -->|
 * 150°C  __________  350°C
 *       |A        B|
 *       |          |
 *       |    0°C   |
 *       |          |
 *  50°C |D________C| 500°C
 */
#define LENGTH 1            // side length, m
#define INIT_TEMP_BODY 0    // °C
#define INIT_TEMP_A 150     // °C
#define INIT_TEMP_B 350     // °C
#define INIT_TEMP_C 500     // °C
#define INIT_TEMP_D 50      // °C



#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define RET_OK        0
#define RET_NEARGS    1
#define RET_BADARGS   2
#define RET_BADMALLOC 3
#define RET_NOTINIT   4

#define NUMARGS       4
#define NPTS_ARG      1
#define NT_ARG        2
#define DT_ARG        3
#define ALPHA_ARG     4

// number of worker threads calculating interior nodes
#define NUM_THREADS   7

#define ERR(err_string)     printf("\e[31mError:\e[0m %s\n", (err_string))
#define WARN(warn_string)   printf("\e[33mWarning:\e[0m %s\n", (warn_string))

#define USAGE printf("\n\e[32mUsage:\e[0m %s [npts] [nt] [dt] [alpha]\n%s\n%s\n%s\n%s\n\n", argv[0],  \
              "  - npts (int): number of grid points in x and y (decimals will be truncated)",        \
              "  - nt (int): number of time steps (decimals will be truncated)",                      \
              "  - dt (float): size of the time steps, in seconds",                                   \
              "  - alpha (float): thermal diffusivity, in m^2/s")

int retval = RET_OK;
#define TEST_RETVAL(value) if((value) != RET_OK){return (value);}

int init_storage(float** a, const int length) {
  *a = (float*)malloc(length * sizeof(float));
  if( *a == NULL ) {
    ERR("could not allocate memory for nodal meshes - try again, or with a smaller value for npts");
    return RET_BADMALLOC;
  } else {
    return RET_OK;
  }
}

// array flipping is 'atomic' (in that I can do it in one line in the main program)
void flip_arrays( float** label_a, float** label_b) {
  float* placeholder = *label_a;
  *label_a = *label_b;
  *label_b = placeholder;
  return;
}

// format: [node 0,0],[node 1,0],...;[node 0,1],[node 1,1],...;...,[node (num_points - 1),(num_points - 1)];\n
int write_data(const int frame_num, const float* array, const int num_points) {
  //printf("%d:", frame_num);
  (void) frame_num;
  for( int y = 0; y < num_points; y++ ) {
    printf("%.2f", array[y * num_points]);
    for( int x = 1; x < num_points; x++ ) {
      printf(",%.2f", array[x + (y * num_points)]);
    }
    printf(";");
  }
  printf("\n");
  return RET_OK;
}

// bitmasks for bool var
#define COPY_MASK     0x01

// for spinning node calculations into other threads
typedef struct {
  unsigned char bools;
  pthread_mutex_t copy_mutex;
  pthread_cond_t copy_sig;
  int row_start;
  int row_end;
  float* current_temps;
  float* previous_temps;
  int npts;
  float fourier;
} calc_data_t;

void* calc_interior( void* calc_data ) {
  calc_data_t* cdp = (calc_data_t*)calc_data;
  pthread_mutex_lock(&cdp->copy_mutex);
  int start = cdp->row_start;
  int end = cdp->row_end;
  cdp->bools |= COPY_MASK;
  pthread_cond_signal(&cdp->copy_sig);
  pthread_mutex_unlock(&cdp->copy_mutex);

  int P;
  for( int y = start; y < end; y++ ) {
    for( int x = 1; x < (cdp->npts - 1); x++) {
      P = x + (y * cdp->npts);
      (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
      (cdp->current_temps)[P] += cdp->fourier * ((cdp->previous_temps)[P + 1] + \
          (cdp->previous_temps)[P - 1] + (cdp->previous_temps)[P + cdp->npts] + \
          (cdp->previous_temps)[P - cdp->npts]);
    }
  }
  return (NULL);
}

void* calc_edges( void* calc_data ) {
  calc_data_t* cdp = (calc_data_t*)calc_data;
  pthread_mutex_lock(&cdp->copy_mutex);
  cdp->bools |= COPY_MASK;
  pthread_cond_signal(&cdp->copy_sig);
  pthread_mutex_unlock(&cdp->copy_mutex);

  int P;
  for( int i = 1; i < (cdp->npts - 1); i++ ) {
    P = i;
    (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
    (cdp->current_temps)[P] += cdp->fourier * ((2 * (cdp->previous_temps)[P + cdp->npts]) + \
        (cdp->previous_temps)[P - 1] + (cdp->previous_temps)[P + 1]);
    P = ((i + 1) * cdp->npts) - 1;
    (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
    (cdp->current_temps)[P] += cdp->fourier * ((2 * (cdp->previous_temps)[P - 1]) +    \
        (cdp->previous_temps)[P - cdp->npts] + (cdp->previous_temps)[P + cdp->npts]);
    P = (cdp->npts * (cdp->npts - 1)) + i;
    (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
    (cdp->current_temps)[P] += cdp->fourier * ((2 * (cdp->previous_temps)[P - cdp->npts]) + \
        (cdp->previous_temps)[P - 1] + (cdp->previous_temps)[P + 1]);
    P = i * cdp->npts;
    (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
    (cdp->current_temps)[P] += cdp->fourier * ((2 * (cdp->previous_temps)[P + 1]) +    \
        (cdp->previous_temps)[P - cdp->npts] + (cdp->previous_temps)[P + cdp->npts]);
  }

  return NULL;
}

void* calc_corners( void* calc_data ) {
  calc_data_t* cdp = (calc_data_t*)calc_data;
  pthread_mutex_lock(&cdp->copy_mutex);
  cdp->bools |= COPY_MASK;
  pthread_cond_signal(&cdp->copy_sig);
  pthread_mutex_unlock(&cdp->copy_mutex);

  int P = 0;
  (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
  (cdp->current_temps)[P] += 2 * cdp->fourier * ((cdp->previous_temps)[P + 1] + (cdp->previous_temps)[P + cdp->npts]);
  P = cdp->npts - 1;
  (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
  (cdp->current_temps)[P] += 2 * cdp->fourier * ((cdp->previous_temps)[P - 1] + (cdp->previous_temps)[P + cdp->npts]);
  P = (cdp->npts * cdp->npts) - 1;
  (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
  (cdp->current_temps)[P] += 2 * cdp->fourier * ((cdp->previous_temps)[P - 1] + (cdp->previous_temps)[P - cdp->npts]);
  P = cdp->npts * (cdp->npts - 1);
  (cdp->current_temps)[P] = (cdp->previous_temps)[P] * (1 - (4 * cdp->fourier));
  (cdp->current_temps)[P] += 2 * cdp->fourier * ((cdp->previous_temps)[P + 1] + (cdp->previous_temps)[P - cdp->npts]);

  return NULL;
}



int main(int argc, char* argv[]) {
  if( argc != (NUMARGS + 1)) {
    ERR("wrong number of args");
    USAGE;
    return RET_NEARGS;
  }

  calc_data_t calc_data;

  calc_data.npts = atoi(argv[NPTS_ARG]);
  const int area = calc_data.npts * calc_data.npts;
  const int nt = atoi(argv[NT_ARG]);
  const float dt = atof(argv[DT_ARG]);
  const float alpha = atof(argv[ALPHA_ARG]);
  // these are all physical values and cannot be zero, so that's a safe check for bad input
  if( (calc_data.npts * nt * dt * alpha) == 0 ) {
    ERR("invalid input - only nonzero numbers are valid");
    USAGE;
    return RET_BADARGS;
  }

  const float dx = (float)LENGTH / calc_data.npts;
  calc_data.fourier = (float)(alpha * dt) / (float)(dx * dx);

  // We need 2 arrays to store the value of each node at the current and previous time step.
  // We'll write the results out to a file as they're calculated so we don't clobber memory.
  // We'll also swap which one is current and previous each iteration so we don't have to copy
  // memory every time. This is probably premature optimization, but it's a concept I've wanted to
  // try out for a while now, so I'm just gonna go for it.
  retval = init_storage(&calc_data.current_temps, area); TEST_RETVAL(retval);
  retval = init_storage(&calc_data.previous_temps, area); TEST_RETVAL(retval);

  // Apply initial conditions - the project assignment seems to imply that the outermost nodes
  // suddenly change from the body's initial temperature to the edge initial temperatures at t=0.
  // Effectively, this simulates a body with a side length of (calc_data.npts - 2)  being suddenly surrounded
  // by material with a thickness of 1 node that has been pre-heated to the specified conditions
  // (neglecting contact resistance). We'll assume the exterior of this material is adiabatic.
  // memcpy is only good for 1 byte data; we need to do it the long way
  for( int y = 0; y < calc_data.npts; y++ ) {
    for( int x = 0; x < calc_data.npts; x++ ) {
      calc_data.current_temps[x + (y * calc_data.npts)] = INIT_TEMP_BODY;
    }
  }
  for( int i = 0; i < calc_data.npts; i++ ) {
    // top
    calc_data.current_temps[i] = INIT_TEMP_A + \
        (i * ((float)(INIT_TEMP_B - INIT_TEMP_A) / (float)(calc_data.npts - 1)));
    // bottom
    calc_data.current_temps[i + ((calc_data.npts - 1) * calc_data.npts)] = \
        INIT_TEMP_D + (i * ((float)(INIT_TEMP_C - INIT_TEMP_D) / (float)(calc_data.npts - 1)));
    // left
    calc_data.current_temps[i * calc_data.npts] = \
        INIT_TEMP_A + (i * ((float)(INIT_TEMP_D - INIT_TEMP_A) / (float)(calc_data.npts - 1)));
    // right
    calc_data.current_temps[((i + 1) * calc_data.npts) - 1] = \
      INIT_TEMP_B + (i * ((float)(INIT_TEMP_C - INIT_TEMP_B) / (float)(calc_data.npts - 1)));
  }

  // calculate - we're assuming an adiabatic exterior
  // worker thread pool
  pthread_t thread_pool[NUM_THREADS];
  // initialize the condtion var to let the main thread know the child has copied the data...
  pthread_cond_init(&calc_data.copy_sig, NULL);
  // ...and the mutex it needs to go with it
  pthread_mutex_init(&calc_data.copy_mutex, NULL);

  for( int i = 0; i < nt; i++ ) {
    write_data(i, calc_data.current_temps, calc_data.npts);
    // there's probably a cleaner way to get the arrays to flip, but that's a later problem
    flip_arrays(&calc_data.current_temps, &calc_data.previous_temps);

    // kick off interior node calc threads
    // Lock to write to data struct, then unlock and wait for the newly created thread to copy
    // in the data. Repeat. After the last thread has copied, we don't need the lock anymore
    // so we have to manually unlock it.
    pthread_mutex_lock(&calc_data.copy_mutex);
    for( int j = 0; j < (NUM_THREADS - 1); j++ ) {
      calc_data.row_start = (j * (calc_data.npts / NUM_THREADS)) + 1;
      calc_data.row_end   = ((j + 1) * (calc_data.npts / NUM_THREADS)) + 1;
      calc_data.bools &= ~COPY_MASK;
      pthread_create(&(thread_pool[j]), NULL, calc_interior, (void*)(&calc_data));
      while( !(calc_data.bools & COPY_MASK) ) {
        pthread_cond_wait(&calc_data.copy_sig, &calc_data.copy_mutex);
      }
    }
    calc_data.row_start = ((NUM_THREADS - 1) * (calc_data.npts / NUM_THREADS)) + 1;
    calc_data.row_end = calc_data.npts - 1;
    calc_data.bools &= ~COPY_MASK;
    pthread_create(&(thread_pool[NUM_THREADS - 1]), NULL, calc_interior, (void*)(&calc_data));
    while( !(calc_data.bools & COPY_MASK) ) {
      pthread_cond_wait(&calc_data.copy_sig, &calc_data.copy_mutex);
    }
    pthread_mutex_unlock(&calc_data.copy_mutex);

    // we'll calculate edges and corners in the main thread, since those don't take long
    calc_edges((void*)(&calc_data));
    calc_corners((void*)(&calc_data));

    // wait for interior node calcs
    for( int j = 0; j < NUM_THREADS; j++ ) {
      pthread_join(thread_pool[j], NULL);
    }
  }
  // put out our last round of calculations
  write_data(nt, calc_data.current_temps, calc_data.npts);

  free(calc_data.current_temps);
  free(calc_data.previous_temps);
  return RET_OK;
}
