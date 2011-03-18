#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

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
    
    int si = sizeof(int);
    int sd = sizeof(double);
    int ind = -si;
    int blens[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint indices[] = {ind += si, ind += si, ind += si, 
                          ind += sd, ind += sd, ind += sd, 
                          ind += sd, ind += sd, ind += sd, sizeof(particle_t)};
    MPI_Datatype oldtypes[] = { MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, 
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB  };

    MPI_Datatype PARTICLE;
    MPI_Type_struct(10, blens, indices, oldtypes, &PARTICLE);
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
    
    // Define local block
    int first = partition_offsets[rank];
    int last = partition_offsets[rank] + partition_sizes[rank];
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    double size = set_size(n);
    if (rank == 0)
    {
        init_particles(n, particles);
    }

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

#if DEBUG
    double times[5]; // TODO: Remove this.
#endif

    // Create a grid for optimizing the interactions
    int gridSize = (size/cutoff) + 1; // TODO: Rounding errors?
    grid_init(grid, gridSize);
    for (int i = 0; i < n; ++i)
    {
        particles[i].grid_coord = particles[i].prev_grid_coord = grid_coord_flat(grid.size, particles[i].x, particles[i].y);
        grid_add(grid, particles[i]);
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
        for (int i = first; i < last; i++)
        {
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
                        apply_force(particles[i], particles[neighbour->particle_id]);
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
        for (int i = first; i < last; i++)
        {
            particles[i].prev_grid_coord = particles[i].grid_coord;

            move(particles[i]);

            particles[i].grid_coord = grid_coord_flat(grid.size, particles[i].x, particles[i].y);
            
            // Re-add the particle if it has changed grid position
            if (particles[i].prev_grid_coord != particles[i].grid_coord)
            {
                if (! grid_remove(grid, particles[i], particles[i].prev_grid_coord))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%p'. Code must be faulty. Blame source writer.\n", &particles[i]);
                    exit(3);
                    return 3;
                }

                grid_add(grid, particles[i]);
            }
        }

        #if DEBUG
        times[1] += (read_timer() - start);
        start = read_timer();
        #endif

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0)
        {
            printf("\n ------- \n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Synchronize particles by heartbeat distribution
        int previous = (rank == 0 ? n_proc - 1 : rank -1);
        int next = (rank == n_proc - 1 ? 0 : rank + 1);

        MPI_Status * stat = (MPI_Status*) malloc(sizeof(MPI_Status));
        int lastblock = rank;
        // Send our particles to the next thread
        for (int i = 1; i < n_proc; ++i)
        {
            // Calculate what partition we're currently processing
            int block = rank - i;
            if (block < 0)
                block = n_proc + block;

            if (block == rank)
                continue;

            for(int i = 0; i < n_proc; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == i)
                {
                    printf("[%d] Heartbeat: Send %d (%d - %d) to %d, Receive %d (%d - %d) from %d\n", rank, lastblock, partition_offsets[lastblock], partition_offsets[lastblock]+partition_sizes[lastblock]-1, next, block, partition_offsets[block], partition_offsets[block]+partition_sizes[block]-1, previous);
                    if (i+1 == n_proc)
                        printf("\n");
                    
                    fflush(stdout);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            
            // Heartbeat communication.
            MPI_Sendrecv(particles + partition_offsets[lastblock], partition_sizes[lastblock], 
                         PARTICLE, next, lastblock, 

                         particles + partition_offsets[block], partition_sizes[block], 
                         PARTICLE, previous, block, MPI_COMM_WORLD, stat);
            
            MPI_Barrier(MPI_COMM_WORLD);
            for(int i = 0; i < n_proc; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == i)
                {
                    int count;
                    MPI_Get_count(stat, PARTICLE, &count);
                    printf("[%d]: Received %d particles. Status: [count: %d, cancelled: %d, Source: %d, Tag: %d, Error: %d]\n", 
                                rank, count, stat->count, stat->cancelled, stat->MPI_SOURCE, stat->MPI_TAG, stat->MPI_ERROR);
                    if (i+1 == n_proc)
                        printf("\n");
                    
                    fflush(stdout);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }

            if (stat->MPI_ERROR != MPI_SUCCESS)
            {
                int tmp;
                fprintf(stdout, "[%d] ERROR %d (count: %d, cancelled: %d, Source: %d, Tag: %d, Error: %d)\n! while recieving block %d(%d - %d) from rank %d.\n", 
                        rank, stat->MPI_ERROR, 
                        block, partition_offsets[lastblock], partition_offsets[lastblock]+partition_sizes[lastblock]-1, previous);
                fflush(stdout);
                printf("IDS: ");
                for (int p = partition_offsets[block]; p < partition_offsets[block]+partition_sizes[block]; ++p)
                    printf("%d ", particles[p].id);
                printf("\n");
                fflush(stdout);
                exit(5);
                return 5;
            }

            for (int p = partition_offsets[block]; p < partition_offsets[block]+partition_sizes[block]; ++p)
            {
                // Re-add the particle if it has changed grid position
                if (particles[p].prev_grid_coord != particles[p].grid_coord)
                {
                    if (! grid_remove(grid, particles[p], particles[p].prev_grid_coord))
                    {
                        fprintf(stdout, "Error: Failed to remove particle '%p'. Code must be faulty. Blame source writer.\n", &particles[p]);
                        exit(4);
                        return 4;
                    }
                    grid_add(grid, particles[p]);
                }
            }
            lastblock = block;
        }
        free(stat);
        MPI_Barrier(MPI_COMM_WORLD);
        
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
    free(particles);
    if (fsave)
        fclose(fsave);
    
    MPI_Finalize();
    
    return 0;
}
/* */