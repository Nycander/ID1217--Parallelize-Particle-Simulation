#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "grid.h"

#define DEBUG 1

//
//  Benchmarking program
//
int main(int argc, char **argv)
{    
    //
    //  Process command line parameters
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
    //  Set up MPI
    //
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //
    //  Allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    grid_t grid;
    
    //
    // Define particle_t in mpi
    //
    int si = sizeof(int);
    int sd = sizeof(double);
    int ind = -si;
    int blens[] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint indices[] = {ind += si, 
                          ind += sd, ind += sd, ind += sd, 
                          ind += sd, ind += sd, ind += sd, sizeof(particle_t)};
    MPI_Datatype oldtypes[] = { MPI_INTEGER, 
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB  };

    MPI_Datatype PARTICLE;
    MPI_Type_struct(8, blens, indices, oldtypes, &PARTICLE);
    MPI_Type_commit(&PARTICLE);
    
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
        grid_add(grid, particles[i]);
    }
    
    //  
    //  Set up the data partitioning across processors
    //
    int rows_per_proc = (gridSize + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc((n_proc+1) * sizeof(int));
    for (int i = 0; i < n_proc+1; i++)
    {
        partition_offsets[i] = Min(i * rows_per_proc, gridSize);
    }
    int *partition_sizes = (int*) malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
    {
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    }
    
    // Define local block
    int first = partition_offsets[rank];
    int last = partition_offsets[rank] + Max(0, partition_sizes[rank] - 1);
    if (rank == 0)
        printf("Grid size = %d\n", grid.size);
    printf("[%d]: Rows: %d -> %d\n", rank, first, last);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

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
            save(fsave, n, particles); // Todo: this will BREAK
        
        #if DEBUG
        double start = read_timer();
        #endif

        //
        //  compute all forces
        //
        int locals[n];
        int local_size = 0;
        // Loop over our part of the grid matrix
        for (int r = first; r <= last; r++)
            for (int c = 0; c < grid.size; ++c)
                // Loop over all particles in this cell
                for(linkedlist_t * current = grid.grid[r * grid.size + c]; 
                    current != 0; 
                    current = current->next)
                {
                    int i = current->particle_id;
                    locals[local_size++] = i;

                    particles[i].ax = particles[i].ay = 0;
                    // Use the grid to traverse neighbours
                    int gx = grid_coord(particles[i].x);
                    int gy = grid_coord(particles[i].y);

                    for(int nx = Max(gx - 1, 0); nx <= Min(gx + 1, grid.size-1); nx++)
                        for(int ny = Max(gy - 1, 0); ny <= Min(gy + 1, grid.size-1); ny++)
                            for(linkedlist_t * neighbour = grid.grid[nx * grid.size + ny]; neighbour != 0; neighbour = neighbour->next)
                                apply_force(particles[i], particles[neighbour->particle_id]);
                }

        
        #if DEBUG
        times[0] += (read_timer() - start);
        start = read_timer();
        #endif

        //printf("[%d] Size before clear: %d\n", rank, grid_size(grid));fflush(stdout);
        // Clear rows first-1 and last+1
        if (first > 0)
        {
            grid_clear_row(grid, first-1);
        }
        if (last < grid.size-1)
        {
            grid_clear_row(grid, last+1);
        }

        //printf("[%d] Size before move: %d\n", rank, grid_size(grid));fflush(stdout);

        // Keep track of what processors we will communicate with
        int sent_to[n_proc];
        memset(sent_to, 0, sizeof(int)*n_proc);

        //
        //  Move particles
        //
        for (int i = 0; i < local_size; ++i)
        {
            int id = locals[i];
            int gc = grid_coord_flat(grid.size, particles[id].x, particles[id].y);

            move(particles[id]);

            int gx = grid_coord(particles[id].x);
            
            // The particle has moved to a neighbour
            if (gx < first || gx > last)
            {
                // Not our problem anymore.
                if (! grid_remove(grid, particles[id], gc))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%d'. Code must be faulty. Blame source writer.\n", id);
                    exit(3);
                    return 3;
                }

                grid_add(grid, particles[id]);

                // Send it to neighbour
                int target = 0;
                while(target < n_proc)
                {

                    int f = partition_offsets[target];
                    int l = partition_offsets[target] + Max(0, partition_sizes[target] - 1);

                    if (gx >= f && gx <= l)
                    {
                        //printf("[%d]: Found target for gx=%d which is %d.\n", rank, gx, target);
                        fflush(stdout);
                        break;
                    }
                    
                    target++;
                }

                sent_to[target] = 1;

                //printf("[%d]: Sending %d to %d.\n", rank, id, target);
                fflush(stdout);
                MPI_Request request;
                MPI_Isend(particles+id, 1, PARTICLE, target, target, MPI_COMM_WORLD, &request);
                continue;
            }
            // The particle has moved to a new internal cell
            if (gc != grid_coord_flat(grid.size, particles[id].x, particles[id].y))
            {
                if (! grid_remove(grid, particles[id], gc))
                {
                    fprintf(stdout, "Error: Failed to remove particle '%d'. Code must be faulty. Blame source writer.\n", id);
                    exit(4);
                    return 4;
                }
                grid_add(grid, particles[id]);
            }

        }
        #if DEBUG
        times[1] += (read_timer() - start);
        start = read_timer();
        #endif

        MPI_Request request;
        // Send first row
        if (first > 0)
        {
            sent_to[rank-1] = 1;
            //printf("[%d] Sending first row %d.\n", rank, first);
            for (int c = 0; c < grid.size; ++c)
                for(linkedlist_t * current = grid.grid[first * grid.size + c]; 
                    current != 0; 
                    current = current->next)
                {
                    //printf("[%d] \tSending %d\n", rank, current->particle_id);
                    MPI_Isend(particles+current->particle_id, 1, PARTICLE, rank-1, rank-1, MPI_COMM_WORLD, &request);
                }
        }
        // Send last row
        if (last < grid.size-1)
        {
            sent_to[rank+1] = 1;
            //printf("[%d] Sending last row %d.\n", rank, last);
            for (int c = 0; c < grid.size; ++c)
                for(linkedlist_t * current = grid.grid[(last) * grid.size + c]; 
                    current != 0; 
                    current = current->next)
                {
                    //printf("[%d] \tSending %d\n", rank, current->particle_id);
                    MPI_Isend(particles+current->particle_id, 1, PARTICLE, rank+1, rank+1, MPI_COMM_WORLD, &request);
                }
        }
        particle_t end;
        end.id = -1;

        // Tell neighbouring rows that we are done.
        if (first > 0)            MPI_Isend(&end, 1, PARTICLE, rank-1, rank, MPI_COMM_WORLD, &request);
        if (last  < grid.size-1)  MPI_Isend(&end, 1, PARTICLE, rank+1, rank, MPI_COMM_WORLD, &request);

        #if DEBUG
        times[2] += (read_timer() - start);
        start = read_timer();
        #endif

        //printf("[%d] Size before receive: %d\n", rank, grid_size(grid));fflush(stdout);

        // Receive new particles from neighbouring threads.
        particle_t new_particle;
        int breakcount = 0;

        // Here be dragons.
        int sent_tos[n_proc];
        MPI_Allreduce(sent_to, sent_tos, n_proc, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
        int breaklimit = sent_tos[rank];

        #if DEBUG
        times[3] += (read_timer() - start);
        start = read_timer();
        #endif
        //printf("[%d] \tGonna wait for %d breaks.\n", rank, breaklimit);
        while(breakcount < breaklimit)
        {
            MPI_Status stat;
            MPI_Recv(&new_particle, 1, PARTICLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

            //printf("[%d] \t\tReceived %d.\n", rank, new_particle.id);
            if (new_particle.id == -1)
            {
                breakcount++;
                continue;
            }

            particles[new_particle.id] = new_particle;

            int gx = grid_coord(new_particle.x);
            if (gx < first-1 || gx > last+1)
            {
                //printf("[%d]: %d < %d || %d > %d .... VA I HELVETE?!?!?!? ELÄNDES SATFLÄSK!  <<<--------- FFFFffffffffuuuuuu\n"
                //        ,rank, gx, first-1, gx, last+1);
                fflush(stdout);
                exit(101);
                return 101;
            }

            grid_add(grid, particles[new_particle.id]);
        }
        //printf("[%d] Size after receive: %d  <----\n", rank, grid_size(grid));fflush(stdout);
        
        #if DEBUG
        times[4] += (read_timer() - start);
        start = read_timer();
        #endif
    }
    simulation_time = read_timer() - simulation_time;
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        printf("n = %d, n_procs = %d, simulation time = %f seconds\n", n, n_proc, simulation_time);

    #if DEBUG
    double result[5];
    MPI_Reduce(times, result, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Force: %f\n", rank, result[0]/n_proc);
        printf("Move: %f\n", rank, result[1]/n_proc);
        printf("Send: %f\n", rank, result[2]/n_proc);
        printf("Scatter: %f\n", rank, result[3]/n_proc);
        printf("Receive: %f\n", rank, result[4]/n_proc);
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