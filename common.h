#ifndef __CS_COMMON_H__
#define __CS_COMMON_H__

#include <vector>

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

#if defined(_WIN32) || defined(_WIN64)
#define fmax max
#define fmin min
#pragma warning (disable:4996)
#define snprintf sprintf_s
#endif

inline int Min( int a, int b ) { return a < b ? a : b; }
inline int Max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;


//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;


typedef struct
{
	int size;
	std::vector<std::vector<int>> v;
} grid_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
double set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );


//
// grid routines
//
void grid_init(grid_t & grid, int gridsize);
void grid_clear(grid_t * grid);
void grid_add(grid_t * grid, particle_t * particle, int pid);
void grid_remove(grid_t * grid, particle_t * particles, int pid);
int grid_coord(double c);

void grid_omp_clear(grid_t * grid);
void grid_omp_populate(grid_t * grid, particle_t * particles, int n);

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
