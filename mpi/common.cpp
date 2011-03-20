#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
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
	long seed = (long)time(NULL);
#ifdef _WIN32
	srand(seed);
#else
	srand48(seed);
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
		long seed = (long)time(NULL);
#ifdef _WIN32
		srand(seed);
		int j = rand()%(n-i);
#else
		srand48(seed);
		int j = lrand48()%(n-i);
#endif
		int k = shuffle[j];
		shuffle[j] = shuffle[n-i-1];

		p[i].id = i;
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
void save(FILE * f, int rank, int n, particle_t *p, int * locals, int local_size, MPI_Datatype PARTICLE)
{
	if (rank == 0)
	{
		static bool first = true;
		if( first )
		{
			fprintf( f, "%d %g\n", n, size );
			first = false;
		}
		// Receive all particles from all threads
		int n_others = n - local_size;
		for (int i = 0; i < n_others; ++i)
		{
			particle_t new_particle;
            MPI_Recv(&new_particle, 1, PARTICLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            p[new_particle.id] = new_particle;
		}
		// Write to file
		for( int i = 0; i < n; i++ )
			fprintf( f, "%g %g\n", p[i].x, p[i].y );
	}
	else
	{
		// Send all particles to rank 0
		for (int i = 0; i < local_size; ++i)
		{
			MPI_Send(p+locals[i], 1, PARTICLE, 0, rank, MPI_COMM_WORLD);
		}
	}
	// We don't want any messages spilling out from this function, so wait for everyone.
    MPI_Barrier(MPI_COMM_WORLD);
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
