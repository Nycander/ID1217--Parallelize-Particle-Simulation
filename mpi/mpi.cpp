#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "common.h"
#include "grid.h"

#define DEBUG 1

//
//  benchmarking program
//
int main(int argc, char **argv)
{    
    //
    //  process command line parameters
    //
    if (find_option(argc, argv, "-h") >= 0)
    {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }
    
    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    grid_t grid;
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);
    
    //  
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc((n_proc+1) * sizeof(int));
    for (int i = 0; i < n_proc+1; i++)
        partition_offsets[i] = Min(i * particle_per_proc, n);
    
    int *partition_sizes = (int*) malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition

    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc(nlocal * sizeof(particle_t));
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    double size = set_size(n);
    if (rank == 0)
    {
        init_particles(n, particles);
    }

    //MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );

#if DEBUG
    double times[5];
#endif

    // Create a grid for optimizing the interactions
    int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
    grid_init(grid, gridSize);
    for (int i = 0; i < n; ++i)
    {
        grid_add(grid, &particles[i]);
    }
   
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++)
    {
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if (fsave && (step%SAVEFREQ) == 0)
            save(fsave, n, particles);
        
        #if DEBUG
        double start = read_timer();
        #endif

        //
        //  compute all forces
        //
        for (int i = 0; i < nlocal; i++)
        {
            local[i].ax = local[i].ay = 0;
            // Use the grid to traverse neighbours
            int gx = grid_coord(local[i].x);
            int gy = grid_coord(local[i].y);

            for(int nx = Max(gx - 1, 0); nx <= Min(gx + 1, grid.size-1); nx++)
            {
                for(int ny = Max(gy - 1, 0); ny <= Min(gy + 1, grid.size-1); ny++)
                {
                    linkedlist_t * neighbour = grid.grid[nx * grid.size + ny];
                    while(neighbour != 0)
                    {
                        apply_force(local[i], *(neighbour->value));
                        neighbour = neighbour->next;
                    }
                }
            }
        }
        
        #if DEBUG
        times[0] += (read_timer() - start);
        start = read_timer();
        #endif

        // Clear grid
        int size = grid.size;
        grid_clear(grid);
        grid_init(grid, size);
        
        #if DEBUG
        times[1] += (read_timer() - start);
        start = read_timer();
        #endif

        //
        //  move particles
        //
        for (int i = 0; i < nlocal; i++)
        {
            move(local[i]);
            grid_add(grid, &local[i]);
        }

        #if DEBUG
        times[2] += (read_timer() - start);
        start = read_timer();
        #endif

        for (int r = n_proc-1; r > 0; --r)
        {
            int target = (rank + r) % n_proc;

            MPI_Request request;
            MPI_Isend(local, nlocal, PARTICLE, target, 0, MPI_COMM_WORLD, &request);
        }

        
        #if DEBUG
        times[3] += (read_timer() - start);
        start = read_timer();
        #endif

        // Receive particles from "the cloud" and add it to the grid
        for (int r = 1; r < n_proc; ++r)
        {
            int target = (rank + r) % n_proc;

            // Receive partion_sizes[r] particles
            particle_t parts[partition_sizes[target]];
           
            #if DEBUG
            double rstart = read_timer();
            #endif

            MPI_Status status; // TODO: Check status?
            MPI_Recv(&parts, partition_sizes[target], PARTICLE, target, 0, MPI_COMM_WORLD, &status);

            #if DEBUG
            printf("\tWaited for proc %d for %f s\n", target, read_timer()-rstart);
            #endif

            for (int i = 0; i < partition_sizes[target]; ++i)
            {
                grid_add(grid, &parts[i]);
            }
        }
        
        #if DEBUG
        times[4] += (read_timer() - start);
        start = read_timer();
        #endif
    }
    simulation_time = read_timer() - simulation_time;
    
    if (rank == 0)
        printf("n = %d, n_procs = %d, simulation time = %f seconds\n", n, n_proc, simulation_time);
    
    #if DEBUG
    double result[5];
    MPI_Reduce(times, result, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Force: %f\n", rank, result[0]/n_proc);
        printf("Clear: %f\n", rank, result[1]/n_proc);
        printf("Move: %f\n", rank, result[2]/n_proc);
        printf("Send: %f\n", rank, result[3]/n_proc);
        printf("Receive: %f\n", rank, result[4]/n_proc);
    }
    #endif

    //
    //  release resources
    //
    free(partition_offsets);
    free(partition_sizes);
    free(local);
    free(particles);
    if (fsave)
        fclose(fsave);
    
    MPI_Finalize();
    
    return 0;
}
