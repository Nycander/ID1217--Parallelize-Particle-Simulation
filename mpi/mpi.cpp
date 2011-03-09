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

        //
        //  move particles
        //
        for (int i = 0; i < nlocal; i++)
        {
            move(local[i]);
        }
        
        #if DEBUG
        times[1] += (read_timer() - start);
        start = read_timer();
        #endif

        int r = (rank + 1) % n_proc;
        MPI_Request request;
        MPI_Isend(local, nlocal, PARTICLE, r, 0, MPI_COMM_WORLD, &request);
        for (int i = 0; i < n_proc -1; ++i)
        {
            // TODO: ISend

            MPI_Recv(&parts[r], partition_sizes[r], PARTICLE, r, 0, MPI_COMM_WORLD, &requests[r]);
            grid_add(grid, &parts[i][j]);

            MPI_Request request;
            MPI_Isend(local, nlocal, PARTICLE, r, 0, MPI_COMM_WORLD, &request);
            int r = (r + 1) % n_proc;
        }
        
        #if DEBUG
        times[2] += (read_timer() - start);
        start = read_timer();
        #endif

        // Clear grid
        int size = grid.size;
        grid_clear(grid);
        grid_init(grid, size);

        
        #if DEBUG
        times[3] += (read_timer() - start);
        start = read_timer();
        #endif

        MPI_Request requests[n_proc];

        particle_t *parts[n_proc];
        // Receive particles from "the cloud" and add it to the grid
        for (int r = 0; r < n_proc; ++r)
        {
            int target = (rank + r) % n_proc;
            MPI_Request tmp;
            requests[target] = tmp;
            // Receive partion_sizes[r] particles
            parts[target] = (particle_t*) malloc(sizeof(particle_t) * partition_sizes[target]);
            printf("time to get suspect stuff... size:%d of target:%d\n", partition_sizes[target], target);
            printf("&parts[target]:%d\n", &parts[target]);
            printf("partition_sizes[target]:%d\n", partition_sizes[target]);
            printf("PARTICLE:%d\n", PARTICLE);
            printf("target:%d\n", target);
            printf("0:%d\n", 0);
            printf("MPI_COMM_WORLD:%d\n", MPI_COMM_WORLD);
            printf("&requests[target]:%d\n", &requests[target]);
            MPI_Irecv(&parts[target], partition_sizes[target], PARTICLE, target, 0, MPI_COMM_WORLD, &requests[target]);
            printf("got suspect stuff.\n");
        }
        int result;
        MPI_Status status;
        for (int read = 0; read < n_proc; ++read)
        {
            int i = 0;
            while(true) //while there is still stuff to do.
            {
                if(requests[i] == NULL) continue;
                MPI_Test(&requests[i], &result, &status);

                if(result) //if we should do stuff with this
                {
                    for(int j = 0; j < partition_sizes[i]; j++)
                    {
                        grid_add(grid, &parts[i][j]);
                    }
                    free(parts[i]);
                    requests[i] = NULL;
                    break;
                }
                i = (i + 1) % n_proc;
            
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
    if (rank == 0)
    {
        printf("Force: %f\n", rank, times[0]);
        printf("Move: %f\n", rank, times[1]);
        printf("Send: %f\n", rank, times[2]);
        printf("Clear: %f\n", rank, times[3]);
        printf("Receive: %f\n", rank, times[4]);
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
