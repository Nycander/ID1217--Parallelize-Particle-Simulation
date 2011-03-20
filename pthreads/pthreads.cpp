#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>

#include "common.h"
#include "grid.h"

#define DEBUG 1

//
//  global variables
int n, nsteps, savefreq, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;
grid_t grid;

#if DEBUG
double tforces = 0.0;
double tmove = 0.0;
#endif

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

pthread_mutex_t task_lock = PTHREAD_MUTEX_INITIALIZER;

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first = Min(  thread_id    * particles_per_thread, n );
    int last  = Min( (thread_id+1) * particles_per_thread, n );

    //printf("Thread %d running, particles %d -> %d.\n", thread_id, first, last);

    // Simulate a number of time steps
    for( int step = 0; step < nsteps; step++ )
    {
        #if DEBUG
        double start = read_timer();
        #endif

        // Compute forces
        for(int i = first; i < last; ++i) 
        {
            // Reset acceleration
            particles[i].ax = particles[i].ay = 0;

            // Use the grid to traverse neighbours
            int gx = grid_coord(particles[i].x);
            int gy = grid_coord(particles[i].y);

            for(int nx = Max(gx - 1, 0); nx <= Min(gx + 1, grid.size-1); nx++)
            {
                for(int ny = Max(gy - 1, 0); ny <= Min(gy + 1, grid.size-1); ny++)
                {
                    linkedlist_t * neighbour = grid.grid[nx * grid.size + ny];
                    while(neighbour != 0)
                    {
                        apply_force(particles[i], *(neighbour->value));
                        neighbour = neighbour->next;
                    }
                }
            }
        }

        pthread_barrier_wait( &barrier );
        #if DEBUG
        if (thread_id == 0)
            tforces += (read_timer() - start);
        start = read_timer();
        #endif

        //  Move particles
        for(int i = first; i < last; i++ ) 
        {
            int gc = grid_coord_flat(grid.size, particles[i].x, particles[i].y);

            move(particles[i]);

            // Re-add the particle if it has changed grid position
            if (gc != grid_coord_flat(grid.size, particles[i].x, particles[i].y))
            {
                if (! grid_remove(grid, &particles[i], gc))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%p'. Code must be faulty. Blame source writer.\n", &particles[i]);
                    exit(3);
                }
                grid_add(grid, &particles[i]);
            }
        }

        pthread_barrier_wait( &barrier );
        #if DEBUG
        if (thread_id == 0)
            tmove += (read_timer() - start);
        #endif

        // Save if necessary
        if( thread_id == 0 && fsave && (step%savefreq) == 0 )
            save( fsave, n, particles );
    }

    #if DEBUG
    if (thread_id == 0)
    {
        // Double-check size of n against grid size
        printf("n = %d, grid.n = %d\n", n, grid_size(grid));
    }
    #endif
    
    return NULL;
}

//
//  benchmarking program
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <int> to set the number of steps in the simulation\n" );
        printf( "-f <int> to set the frequency of saving particle coordinates (e.g. each ten's step)\n" );
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    nsteps = read_int(argc, argv, "-s", NSTEPS);
    savefreq = read_int(argc, argv, "-f", SAVEFREQ);
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    

    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double size = set_size( n );
    init_particles( n, particles );

    // Create a grid for optimizing the interactions
    int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
    grid_init(grid, gridSize);
    for (int i = 0; i < n; ++i)
    {
        grid_add(grid, &particles[i]);
    }

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );


    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    #if DEBUG
    printf("Reading: %f\n", tforces);
    printf("Writing: %f\n", tmove);
    #endif

    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    grid_clear(grid);
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
