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

int init_storage(float** a, float** b, const int length) {
  *a = (float*)malloc(length * sizeof(float));
  *b = (float*)malloc(length * sizeof(float));
  if( *a == NULL || *b == NULL ) {
    ERR("could not allocate memory for nodal meshes - try again, or with a smaller value for npts");
    return RET_BADMALLOC;
  } else {
    return RET_OK;
  }
}

// I want to make array swapping 'atomic', in the sense that I don't have to manually keep track of
// which label refers to which data array. This (stateless!) function removes complication - I can
// use a single line to flip the arrays.
void flip_arrays( float** label_a, float** label_b) {
  float* placeholder = *label_a;
  *label_a = *label_b;
  *label_b = placeholder;
  return;
}

// format: [node 0,0],[node 1,0],...;[node 0,1],[node 1,1],...;...,[node (num_points - 1),(num_points - 1)];\n
// add frame num back in if we need it; as of now, its just getting in the way
// data density is very low here, but storage is cheap and readability is more important - I'll
// probably use Python to visualize the data
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

// for spinning interior node calculations into other threads
struct calc_data_struct {
  pthread_cond_t* copy_sig;
  pthread_mutex_t* copy_mutex;
  // for debugging
  int print;
  int row_start;
  int row_end;
  float* current_temps;
  float* previous_temps;
  int npts;
  float fourier;
};
typedef struct calc_data_struct calc_data_t;

void* calc_interior( void* calc_data ) {
  pthread_mutex_lock(((calc_data_t*)calc_data)->copy_mutex);
  calc_data_t data = *((calc_data_t *)calc_data);
  pthread_mutex_unlock(((calc_data_t*)calc_data)->copy_mutex);
  pthread_cond_signal(((calc_data_t*)calc_data)->copy_sig);
  if( data.print ) {
    fprintf(stderr, "[chld] s: %02d, e: %02d, c: %p, p: %p, n: %d, f: %f\n",
        data.row_start, data.row_end, (void*)(data.current_temps),
        (void*)(data.previous_temps), data.npts, data.fourier);
  }
  int P;
  for( int y = data.row_start; y < data.row_end; y++ ) {
    for( int x = 1; x < (data.npts - 1); x++) {
      P = x + (y * data.npts);
      (data.current_temps)[P] = (data.previous_temps)[P] * (1 - (4 * data.fourier));
      (data.current_temps)[P] += data.fourier * ((data.previous_temps)[P + 1] + \
          (data.previous_temps)[P - 1] + (data.previous_temps)[P + data.npts] + \
          (data.previous_temps)[P - data.npts]);
    }
  }
  return (NULL);
}

int main(int argc, char* argv[]) {
  if( argc != (NUMARGS + 1)) {
    ERR("wrong number of args");
    USAGE;
    return RET_NEARGS;
  }

  const int npts = atoi(argv[NPTS_ARG]);
  const int nt = atoi(argv[NT_ARG]);
  const float dt = atof(argv[DT_ARG]);
  const float alpha = atof(argv[ALPHA_ARG]);
  // these are all physical values and cannot be zero, so that's a safe check for bad input
  if( (npts * nt * dt * alpha) == 0 ) {
    ERR("invalid input - only nonzero numbers are valid");
    USAGE;
    return RET_BADARGS;
  }

  const float dx = (float)LENGTH / npts;
  const float fourier = (float)(alpha * dt) / (float)(dx * dx);

  // We need 2 arrays to store the value of each node at the current and previous time step.
  // We'll write the results out to a file as they're calculated so we don't clobber memory.
  // We'll also swap which one is current and previous each iteration so we don't have to copy
  // memory every time. This is probably premature optimization, but it's a concept I've wanted to
  // try out for a while now, so I'm just gonna go for it.
  float* arr_a;
  float* arr_b;
  retval = init_storage(&arr_a, &arr_b, npts * npts); TEST_RETVAL(retval);
  float* current_temps = arr_a;
  float* previous_temps = arr_b;

  // Apply initial conditions - the project assignment seems to imply that the outermost nodes
  // suddenly change from the body's initial temperature to the edge initial temperatures at t=0.
  // Effectively, this simulates a body with a side length of (npts - 2)  being suddenly surrounded
  // by material with a thickness of 1 node that has been pre-heated to the specified conditions
  // (neglecting contact resistance). We'll assume the exterior of this material is adiabatic.
  // memcpy is only good for 1 byte data; we need to do it the long way
  for( int y = 0; y < npts; y++ ) {
    for( int x = 0; x < npts; x++ ) {
      current_temps[x + (y * npts)] = INIT_TEMP_BODY;
    }
  }
  for( int i = 0; i < npts; i++ ) {
    // top
    current_temps[i] = INIT_TEMP_A + (i * ((float)(INIT_TEMP_B - INIT_TEMP_A) / (float)(npts - 1)));
    // bottom
    current_temps[i + ((npts - 1) * npts)] = INIT_TEMP_D + (i * ((float)(INIT_TEMP_C - INIT_TEMP_D) / (float)(npts - 1)));
    // left
    current_temps[i * npts] = INIT_TEMP_A + (i * ((float)(INIT_TEMP_D - INIT_TEMP_A) / (float)(npts - 1)));
    // right
    current_temps[((i + 1) * npts) - 1] = INIT_TEMP_B + (i * ((float)(INIT_TEMP_C - INIT_TEMP_B) / (float)(npts - 1)));
  }

  // calculate - we're assuming an adiabatic exterior
  // save space in equations by pre-calculating the point coordinate
  int P;
  // prepare struct to pass data to worker threads
  calc_data_t calc_data;
  calc_data.print = 0;
  calc_data.npts = npts;
  calc_data.fourier = fourier;
  // worker thread pool
  pthread_t thread_pool[NUM_THREADS];
  // initialize the condtion var to let the main thread know the child has copied the data...
  pthread_cond_t data_change_cond;
  calc_data.copy_sig = &data_change_cond;
  pthread_cond_init(&(data_change_cond), NULL);
  // ...and the mutex it needs to go with it
  pthread_mutex_t data_change_mutex;
  calc_data.copy_mutex = &data_change_mutex;
  pthread_mutex_init(&(data_change_mutex), NULL);
  // we need to take this lock outside the loop, explained below
  pthread_mutex_lock(calc_data.copy_mutex);
  for( int i = 0; i < nt; i++ ) {
    write_data(i, current_temps, npts);
    // there's probably a cleaner way to get the arrays to flip, but that's a later problem
    flip_arrays(&current_temps, &previous_temps);
    calc_data.current_temps = current_temps;
    calc_data.previous_temps = previous_temps;

    // kick off interior node calc threads
    // We manually locked the mutex before the loop so we know we can write to the calc_data struct.
    // Once we've done that, we unlock and wait for the most recently created thread to copy its
    // data, then re-lock and update with the next thread's info. Wash, rinse, repeat - after the
    // last thread signals that it's done copying, we lock the mutex and hold that lock until the
    // loop starts again (that's why the initial lock is out of the loop - otherwise we'd deadlock
    // trying to acquire a lock we already have).
    for( int j = 0; j < (NUM_THREADS - 1); j++ ) {
      calc_data.row_start = (j * (npts / NUM_THREADS)) + 1;
      calc_data.row_end   = ((j + 1) * (npts / NUM_THREADS)) + 1;
      pthread_create(&(thread_pool[j]), NULL, calc_interior, (void*)(&calc_data));
      pthread_cond_wait(calc_data.copy_sig, calc_data.copy_mutex);
    }
    calc_data.row_start = ((NUM_THREADS - 1) * (npts / NUM_THREADS)) + 1;
    calc_data.row_end = npts - 1;
    pthread_create(&(thread_pool[NUM_THREADS - 1]), NULL, calc_interior, (void*)(&calc_data));
    pthread_cond_wait(calc_data.copy_sig, calc_data.copy_mutex);

    /*
    fprintf(stderr, "[prnt] s: %02d, e: %02d, c: %p, p: %p, n: %d, f: %f\n",
        calc_data.row_start, calc_data.row_end, (void*)(calc_data.current_temps),
        (void*)(calc_data.previous_temps), calc_data.npts, calc_data.fourier);

    // hardcoding for 2 threads to debug
    calc_data.row_start = 1; calc_data.row_end = 21; calc_data.print = 1;
    pthread_create(&(thread_pool[0]), NULL, calc_interior, (void*)(&calc_data));
    pthread_cond_wait(calc_data.copy_sig, calc_data.copy_mutex);
    calc_data.row_start = 21; calc_data.row_end = 39; calc_data.print = 1;
    pthread_create(&(thread_pool[1]), NULL, calc_interior, (void*)(&calc_data));
    pthread_cond_wait(calc_data.copy_sig, calc_data.copy_mutex);
    */

    // deal with edges
    for( int j = 1; j < (npts - 1); j++ ) {
      P = j;
      current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
      current_temps[P] += fourier * ((2 * previous_temps[P + npts]) + \
          previous_temps[P - 1] + previous_temps[P + 1]);
      P = ((j + 1) * npts) - 1;
      current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
      current_temps[P] += fourier * ((2 * previous_temps[P - 1]) +    \
          previous_temps[P - npts] + previous_temps[P + npts]);
      P = (npts * (npts - 1)) + j;
      current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
      current_temps[P] += fourier * ((2 * previous_temps[P - npts]) + \
          previous_temps[P - 1] + previous_temps[P + 1]);
      P = j * npts;
      current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
      current_temps[P] += fourier * ((2 * previous_temps[P + 1]) +    \
          previous_temps[P - npts] + previous_temps[P + npts]);
    }

    // deal with corners
    P = 0;
    current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
    current_temps[P] += 2 * fourier * (previous_temps[P + 1] + previous_temps[P + npts]);
    P = npts - 1;
    current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
    current_temps[P] += 2 * fourier * (previous_temps[P - 1] + previous_temps[P + npts]);
    P = (npts * npts) - 1;
    current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
    current_temps[P] += 2 * fourier * (previous_temps[P - 1] + previous_temps[P - npts]);
    P = npts * (npts - 1);
    current_temps[P] = previous_temps[P] * (1 - (4 * fourier));
    current_temps[P] += 2 * fourier * (previous_temps[P + 1] + previous_temps[P - npts]);

    // wait for interior node calcs
    for( int j = 0; j < NUM_THREADS; j++ ) {
      pthread_join(thread_pool[j], NULL);
    }
  }
  // put out our last round of calculations
  write_data(nt, current_temps, npts);

  free(arr_a);
  free(arr_b);
  return RET_OK;
}
