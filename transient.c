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

#define ERR(err_string)     printf("\e[31mError:\e[0m %s\n", (err_string))
#define WARN(warn_string)   printf("\e[33mWarning:\e[0m %s\n", (warn_string))

#define USAGE printf("\n\e[32mUsage:\e[0m %s [npts] [nt] [dt] [alpha]\n%s\n%s\n%s\n%s\n\n", argv[0],  \
              "  - npts (int): number of grid points in x and y (decimals will be truncated)",        \
              "  - nt (int): number of time steps (decimals will be truncated)",                      \
              "  - dt (float): size of the time steps, in seconds",                                   \
              "  - alpha (float): thermal diffusivity, in m^2/s")

#define TEST_RETVAL(value) if((value) != RET_OK){return (value);}

// for testing return values
int retval = RET_OK;

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

int write_to_file(const float* array, const int num_points, const int fd) {
  return RET_OK;
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

//  /* array flip debug
  for( int y = 0; y < npts; y++ ) {
    for( int x = 0; x < npts; x++ ) {
      printf("%.1f ", current_temps[x + (y * npts)]);
    }
    printf("\n");
  }
  flip_arrays(&current_temps, &previous_temps);
  for( int y = 0; y < npts; y++ ) {
    for( int x = 0; x < npts; x++ ) {
      printf("%.1f ", current_temps[x + (y * npts)]);
    }
    printf("\n");
  }
  flip_arrays(&previous_temps, &current_temps);
  for( int y = 0; y < npts; y++ ) {
    for( int x = 0; x < npts; x++ ) {
      printf("%.1f ", current_temps[x + (y * npts)]);
    }
    printf("\n");
  }
//  */

  for( int i = 0; i < nt; i++ ) {
    write_to_file(current_temps, npts, 0);
    flip_arrays(&current_temps, &previous_temps);
    // calculate
  }

  free(arr_a);
  free(arr_b);
  return RET_OK;
}
