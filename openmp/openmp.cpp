#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "common.h"
#include "grid.h"

#define DEBUG 1

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
    int n_threads = read_int(argc, argv, "-p", 2);
    char *savename = read_string( argc, argv, "-o", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double size = set_size( n );
    init_particles( n, particles );

    // Create a grid for optimizing the interactions
    int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
    grid_t grid;
    grid_init(grid, gridSize);
    for (int i = 0; i < n; ++i)
    {
        grid_add(grid, &particles[i]);
    }

    // Set number of threads
    omp_set_num_threads(n_threads);

    // Simulate a number of time steps
    double simulation_time = read_timer( );

    double times[5];
    for (int i = 0; i < 5; ++i)
    {
        times[i] = 0.0;
    }

    #pragma omp parallel
    for( int step = 0; step < NSTEPS; step++ )
    {
        #if DEBUG
        double start;
        #pragma omp single
        {
            start = omp_get_wtime();
        }
        #endif

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
                    linkedlist_t * curr = grid.grid[x * grid.size + y];
                    while(curr != 0)
                    {
                        apply_force(particles[i], *(curr->value));
                        curr = curr->next;
                    }
                }
            }
        }

        #if DEBUG
        #pragma omp single
        {
            times[0] += (omp_get_wtime() - start);
            start = omp_get_wtime();
        }
        #endif
        
        // Clear grid
        #pragma omp single
        {
            int size = grid.size;
            grid_clear(grid);
            grid_init(grid, size);
        }

        #if DEBUG
        #pragma omp single
        {
            times[1] += (omp_get_wtime() - start);
            start = omp_get_wtime();
        }
        #endif

        // Move particles
        #pragma omp for
        for(int i = 0; i < n; i++) 
            move(particles[i]);

        #if DEBUG
        #pragma omp single
        {
            times[2] += (omp_get_wtime() - start);
            start = omp_get_wtime();
        }
        #endif

        // Re-populate grid
        #pragma omp single
        for(int i = 0; i < n; i++)
            grid_add(grid, &particles[i]);
        
        #if DEBUG
        #pragma omp single
        {
            times[3] += (omp_get_wtime() - start);
            start = omp_get_wtime();
        }
        #endif

        // Save if necessary
        #pragma omp master
        if(fsave && (step%SAVEFREQ) == 0)
            save(fsave, n, particles);

        #if DEBUG
        #pragma omp single
        {
            times[4] += (omp_get_wtime() - start);
            start = omp_get_wtime();
        }
        #endif
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf("n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );

    #if DEBUG
    printf("Forces: %f s\n", times[0]);
    printf("Clear:  %f s\n", times[1]);
    printf("Moving: %f s\n", times[2]);
    printf("Re-pop: %f s\n", times[3]);
    printf("Saving: %f s\n", times[4]);
    #endif

    grid_clear(grid);
    free(particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
