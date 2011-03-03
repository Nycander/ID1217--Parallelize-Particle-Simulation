#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double size = set_size( n );
    init_particles( n, particles );

    // Create a grid for optimizing the interactions
    int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
    std::vector<int> tmp[gridSize*gridSize];
    grid_t grid;
    grid.v = tmp;
    grid_init(&grid, gridSize);

    // Simulate a number of time steps
    double simulation_time = read_timer( );

    #pragma omp parallel

    #pragma omp for
    for (int i = 0; i < n; ++i)
    {
        grid_add(&grid, &particles[i], i);
    }

    for( int step = 0; step < 1000; step++ )
    {
        // Compute all forces
        #pragma omp for
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
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );

        // Reset grid
        grid_omp_clear(&grid);

        // Re-populate grid
        #pragma omp for
        for(int i = 0; i < n; i++)
            grid_add(&grid, &particles[i], i);
        
        // Save if necessary
        #pragma omp master
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf("n = %d,\tsimulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
