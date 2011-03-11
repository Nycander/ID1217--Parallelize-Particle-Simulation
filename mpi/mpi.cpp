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
    MPI_Type_contiguous(1, MPI_INTEGER, &PARTICLE);
    MPI_Type_commit(&PARTICLE);
    
    //  
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc((n_proc+1) * sizeof(int));

    for (int i = 0; i < n_proc+1; i++)
    {
        partition_offsets[i] = Min(i * particle_per_proc, n);
    }
    
    int *partition_sizes = (int*) malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
    {
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    }
    
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
        particles[i].prevgrid = grid_coord_flat(grid.size, particles[i].x, particles[i].y);
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

        //
        //  move particles
        //
        for (int i = 0; i < nlocal; i++)
        {
            int gc = grid_coord_flat(grid.size, local[i].x, local[i].y);

            local[i].prevgrid = gc;

            move(local[i]);

            // Re-add the particle if it has changed grid position
            if (gc != grid_coord_flat(grid.size, local[i].x, local[i].y))
            {
                if (! grid_remove(grid, particles + partition_offsets[rank] + i, gc))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%p'. Code must be faulty. Blame source writer.\n", &local[i]);
                    exit(3);
                    return 3;
                }

                grid_add(grid, particles + partition_offsets[rank] + i);
            }
        }

        #if DEBUG
        times[1] += (read_timer() - start);
        start = read_timer();
        #endif

        // Synchronize particles by pipeline
        int r = (rank + 1) % n_proc;

        // Get all other particle blocks
        int previous = (rank == 0 ? n_proc - 1 : rank -1);
        int next = (rank == n_proc - 1 ? 0 : rank + 1);


        // Send our particles to the next thread
        MPI_Request request;
        // TODO: Let the 'tag' be the block id.
        if (n_proc > 1) MPI_Isend(local, nlocal, PARTICLE, next, rank, MPI_COMM_WORLD, &request);

        for (int i = 0; i < n_proc; ++i)
        {
            int block = rank - i;
            if (block < 0)
                block = n_proc + block;

            if (block == rank)
                continue;
            
            // TODO: Let the 'tag' be the block id.
            MPI_Status * stat = (MPI_Status*) malloc(sizeof(MPI_Status));
            MPI_Recv(particles + partition_offsets[block], 
                     partition_sizes[block], 
                     PARTICLE, previous, block, MPI_COMM_WORLD, stat);
            
            if (stat->MPI_ERROR != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR! ERROR! %d\n", stat->MPI_ERROR);
                exit(5);
                return 5;
            }

            for (int p = 0; p < partition_sizes[block]; ++p)
            {
                particle_t * part = particles+partition_offsets[block]+p;
                // Re-add the particle if it has changed grid position
                if (part->prevgrid != grid_coord_flat(grid.size, part->x, part->y))
                {
                    printf("%d: Removing %p ....\n", rank, part);

                    printf("Before %d: [", rank);
                    linkedlist_t * curr = grid.grid[part->prevgrid];
                    while(curr != 0)
                    {
                        printf("%p ", curr->value);
                        curr = curr->next;
                    }
                    printf("]\n");
                    fflush(stdout);

                    if (! grid_remove(grid, part, part->prevgrid))
                    {
                        fprintf(stdout, "Error: Failed to remove particle '%p'. Code must be faulty. Blame source writer.\n", &particles[i]);
                        exit(4);
                        return 4;
                    }
                    printf("REMOVED %d COCK(S) FROM ASS\n", rank); fflush(stdout);

                    printf("After %d: [", rank);
                    curr = grid.grid[part->prevgrid];
                    while(curr != 0)
                    {
                        printf("%p ", curr->value);
                        curr = curr->next;
                    }
                    printf("]\n");
                    fflush(stdout);

                    grid_add(grid, part);
                }
            }

            if (block == next)
                continue;

            MPI_Request req;
            MPI_Isend(particles + partition_offsets[block], partition_sizes[block], 
                      PARTICLE, next, block, MPI_COMM_WORLD, &req);
        }
        
        #if DEBUG
        times[2] += (read_timer() - start);
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
