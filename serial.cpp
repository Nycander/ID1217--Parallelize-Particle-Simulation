#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "common.h"
#include <vector>

#define DEBUG 0
using namespace std;

int gridCoord(double c)
{
	return c / cutoff; // TODO: Roundoff errors?
}

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
	
	FILE *fsave = savename ? fopen( savename, "w" ) : stdout ;
	particle_t * particles = (particle_t*) malloc( n * sizeof(particle_t) );
	double size = set_size( n );
	init_particles( n, particles );
	// Create a grid for optimizing the interactions
	int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
	vector<int> grid[gridSize][gridSize];

	for(int i = 0; i < gridSize; i++)
	for(int j = 0; j < gridSize; j++)
		grid[i][j] = vector<int>();

	printf("Creating grid of size %dx%d...\n", gridSize, gridSize); fflush(stdout);

	for(int i = 0; i < n; i++)
	{
		particle_t * p = &particles[i];
		int gridx = gridCoord(p->x);
		int gridy = gridCoord(p->y);
		grid[gridx][gridy].push_back(i);
	}
	

	//
	//  simulate a number of time steps
	//
	double simulation_time = read_timer( );
	for( int step = 0; step < NSTEPS; step++ )
	{
		//
		//  compute forces
		//
		for( int i = 0; i < n; i++ )
		{
			 // Reset acceleration
			particles[i].ax = particles[i].ay = 0;

			int gx = gridCoord(particles[i].x);
			int gy = gridCoord(particles[i].y);
#if DEBUG
			printf("\tForcecheck: %.3f,%.3f (grid: %d,%d)...\n", particles[i].x, particles[i].y, gx, gy); fflush(stdout);
#endif
			for(int x = Max(gx - 1, 0); x <= Min(gx + 1, gridSize-1); x++)
			{
				for(int y = Max(gy - 1, 0); y <= Min(gy + 1, gridSize-1); y++)
				{
                    for(int p = 0; p < grid[x][y].size(); p++)
                    {
                        apply_force(particles[i], particles[grid[x][y][p]]);
                    }
				}
			}
		}

#if DEBUG
		printf("\nMoving particles phase...\n");
		fflush(stdout);
#endif
		// Reset grid
		for(int i = 0; i < gridSize; i++)
		for(int j = 0; j < gridSize; j++)
			grid[i][j].clear();

		//
		//  move particles
		//
		for( int i = 0; i < n; i++ ) 
		{
#if DEBUG
			
			printf("\tMove: %.3f,%.3f (%d,%d) -> ", particles[i].x, particles[i].y, gridCoord(particles[i].x), gridCoord(particles[i].y));
			fflush(stdout);
#endif
			move( particles[i] );

			int gridx = gridCoord(particles[i].x);
			int gridy = gridCoord(particles[i].y);

			grid[gridx][gridy].push_back(i);
#if DEBUG
			printf("%.3f,%.3f (%d,%d)\n", particles[i].x, particles[i].y, gridx, gridy);
			fflush(stdout);
#endif
		}

#if DEBUG
		printf("\nSaving?\n"); fflush(stdout);
#endif		

		//
		//  save if necessary
		//
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
