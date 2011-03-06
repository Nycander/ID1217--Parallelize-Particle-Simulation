#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#ifndef _WIN32
#include <sys/time.h>
#else
#include "gettimeofday.h"
#endif

#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
	static bool initialized = false;
	static struct timeval start;
	struct timeval end;
	if( !initialized )
	{
		gettimeofday( &start, NULL );
		initialized = true;
	}
	gettimeofday( &end, NULL );
	return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
double set_size( int n )
{
	size = sqrt( density * n );
	return size;
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
// TODO: Use portable c code
#ifdef _WIN32
	srand((long)time(NULL));
#else
	srand48( (long)time(NULL) );
#endif

	int sx = (int)ceil(sqrt((double)n));
	int sy = (n+sx-1)/sx;

	int *shuffle = (int*)malloc( n * sizeof(int) );
	for( int i = 0; i < n; i++ )
		shuffle[i] = i;

	for( int i = 0; i < n; i++ ) 
	{
		//
		//  make sure particles are not spatially sorted
		//
#ifdef _WIN32
		srand((long)time(NULL));
		int j = rand()%(n-i);
#else
		srand48((long)time(NULL));
		int j = lrand48()%(n-i);
#endif
		int k = shuffle[j];
		shuffle[j] = shuffle[n-i-1];

		//
		//  distribute particles evenly to ensure proper spacing
		//
		p[i].x = size*(1.+(k%sx))/(1+sx);
		p[i].y = size*(1.+(k/sx))/(1+sy);

		//
		//  assign random velocities within a bound
		//
#ifdef _WIN32
		p[i].vx = (double)(rand()/(double)INT_MAX)*2-1;
		p[i].vy = (double)(rand()/(double)INT_MAX)*2-1;
#else
		p[i].vx = drand48()*2-1;
		p[i].vy = drand48()*2-1;
#endif
	}
	free( shuffle );
}

//
// initialize grid and fill it with particles
// 
void grid_init(grid_t & grid, int size)
{
	grid.size = size;
	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < size; j++)
		{
			std::vector<int> hej;
			grid.v[(i * size) + j] = hej;
		}
	}
}

//
// populates a grid with particles openmp style
//
void grid_omp_populate(grid_t * grid, particle_t * particles, int n)
{
	#pragma omp for
	for (int i = 0; i < n; ++i)
	{
		int gridx = grid_coord(particles[i].x);
		int gridy = grid_coord(particles[i].y);

		#pragma omp critical
		grid->v[gridx * grid->size + gridy].push_back(i);
	}
}

void grid_add(grid_t & grid, particle_t * p, int pid)
{
	int gridx = grid_coord(p->x);
	int gridy = grid_coord(p->y);

	grid.v[gridx * grid.size + gridy].push_back(pid);
}

//
// grid move
//
void grid_remove(grid_t & grid, particle_t * p, int pid)
{
	int grid_x = grid_coord(p->x);
	int grid_y = grid_coord(p->y);
	int gridCoord = grid_x * grid.size + grid_y;

	// Remove particle from grid slot
	particle_t * particle;
	for (unsigned int i = 0; i < grid.v[gridCoord].size(); ++i)
	{
		if (grid.v[gridCoord][i] == pid)
		{
			grid.v[gridCoord].erase(grid.v[gridCoord].begin() + i, grid.v[gridCoord].begin() + i + 1);
			break;
		}
	}
}

//
// clears a grid from values
//
void grid_clear(grid_t * grid)
{
	for(int i = 0; i < grid->size; i++)
		for(int j = 0; j < grid->size; j++)
			grid->v[i * grid->size + j].clear();
}
//
// clears a grid from values openmp style
//
void grid_omp_clear(grid_t * grid)
{
	#pragma omp for
	for(int i = 0; i < grid->size; i++)
	for(int j = 0; j < grid->size; j++)
		grid->v[i * grid->size + j].clear();
}

//
// Calculate the grid coordinate from a real coordinate
//
int grid_coord(double c)
{
	return (int)(c / cutoff);
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{
	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;
	if( r2 > cutoff*cutoff )
		return;
	r2 = fmax( r2, min_r*min_r );
	double r = sqrt( r2 );

	//
	//  very simple short-range repulsive force
	//
	double coef = ( 1 - cutoff / r ) / r2 / mass;
	particle.ax += coef * dx;
	particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
	//
	//  slightly simplified Velocity Verlet integration
	//  conserves energy better than explicit Euler method
	//
	p.vx += p.ax * dt;
	p.vy += p.ay * dt;
	p.x  += p.vx * dt;
	p.y  += p.vy * dt;

	//
	//  bounce from walls
	//
	while( p.x < 0 || p.x > size )
	{
		p.x  = p.x < 0 ? -p.x : 2*size-p.x;
		p.vx = -p.vx;
	}
	while( p.y < 0 || p.y > size )
	{
		p.y  = p.y < 0 ? -p.y : 2*size-p.y;
		p.vy = -p.vy;
	}
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
	static bool first = true;
	if( first )
	{
		fprintf( f, "%d %g\n", n, size );
		first = false;
	}
	for( int i = 0; i < n; i++ )
		fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
	for( int i = 1; i < argc; i++ )
		if( strcmp( argv[i], option ) == 0 )
			return i;
	return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
	int iplace = find_option( argc, argv, option );
	if( iplace >= 0 && iplace < argc-1 )
		return atoi( argv[iplace+1] );
	return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
	int iplace = find_option( argc, argv, option );
	if( iplace >= 0 && iplace < argc-1 )
		return argv[iplace+1];
	return default_value;
}
