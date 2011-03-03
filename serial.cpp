#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <vector>

#include "common.h"

#define DEBUG 0
using namespace std;

//
//  benchmarking program
//
int main( int argc, char **argv )
{
	if( find_option( argc, argv, "-h" ) >= 0 )
	{
		printf( "Options:\n" );
		printf( "-h to see this help\n" );
		printf( "-n <int> to set the number of particles\n" );
		printf( "-o <filename> to specify the output file name\n" );
		return 0;
	}
	
	int n = read_int( argc, argv, "-n", 1000 );

	char *savename = read_string( argc, argv, "-o", NULL );
	
	FILE *fsave = savename ? fopen( savename, "w" ) : NULL ;
	particle_t * particles = (particle_t*) malloc( n * sizeof(particle_t) );
	double size = set_size( n );
	init_particles( n, particles );


	// Create a grid for optimizing the interactions
	int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
	std::vector<int> tmp[gridSize*gridSize];
	grid_t grid;
	grid.v = tmp;
	grid_init(&grid, gridSize);
	for (int i = 0; i < n; ++i)
	{
		grid_add(&grid, &particles[i], i);
	}
	
	// Simulate a number of time steps
	double simulation_time = read_timer( );
	for( int step = 0; step < NSTEPS; step++ )
	{
		// Compute forces
		for( int i = 0; i < n; i++ )
		{
			// Reset acceleration
			particles[i].ax = particles[i].ay = 0;

            // Use the grid to traverse neighbours
			int gx = grid_coord(particles[i].x);
			int gy = grid_coord(particles[i].y);

			for(int x = Max(gx - 1, 0); x <= Min(gx + 1, gridSize-1); x++)
			{
				for(int y = Max(gy - 1, 0); y <= Min(gy + 1, gridSize-1); y++)
				{
                    for(int p = 0; p < grid.v[x * grid.size + y].size(); p++)
                    {
                        apply_force(particles[i], particles[grid.v[x * grid.size + y][p]]);
                    }
				}
			}
		}
	
		// Move particles
		for( int i = 0; i < n; i++ ) 
		{
			move( particles[i] );
		}

		// Reset grid
		grid_clear(&grid);

		// Re-populate grid
		for (int i = 0; i < n; ++i)
		{
			grid_add(&grid, &particles[i], i);
		}

		// Save if necessary
		if( fsave && (step%SAVEFREQ) == 0 )
			save( fsave, n, particles );
	}
	simulation_time = read_timer( ) - simulation_time;
	
	printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
	
	free( particles );
	if( fsave )
		fclose( fsave );
	
	return 0;
}
