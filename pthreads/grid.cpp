#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include "grid.h"


//
// Calculate the grid coordinate from a real coordinate
//
int grid_coord(double c)
{
    return (int)floor(c / cutoff);
}


//
// initialize grid and fill it with particles
// 
void grid_init(grid_t & grid, int size)
{
    grid.size = size;

    // Initialize grid
    grid.grid = (linkedlist**) malloc(sizeof(linkedlist*) * size * size);

    if (grid.grid == NULL)
    {
        fprintf(stderr, "Error: Could not allocate memory for the grid!\n");
        exit(1);
    }

    memset(grid.grid, 0, sizeof(linkedlist*) * size * size);

    // Initialize locks
    grid.lock = (pthread_mutex_t*) malloc(sizeof(pthread_mutex_t) * size * size);

    if (grid.lock == NULL)
    {
        fprintf(stderr, "Error: Could not allocate memory for the locks!\n");
        exit(2);
    }

    for (int i = 0; i < size*size; ++i)
    {
        pthread_mutex_init(&grid.lock[i], NULL);
    }
}

//
// adds a particle pointer to the grid
//
void grid_add(grid_t & grid, particle_t * p)
{
    int gridx = grid_coord(p->x);
    int gridy = grid_coord(p->y);

    int gridCoord = gridx * grid.size + gridy;

    linkedlist_t * newElement = (linkedlist_t *) malloc(sizeof(linkedlist));
    newElement->value = p;

    // Beginning of critical section
    pthread_mutex_lock(&grid.lock[gridCoord]);
    linkedlist_t * tmp = grid.grid[gridCoord];

    newElement->next = tmp;

    grid.grid[gridCoord] = newElement;
    // End of critical section
    pthread_mutex_unlock(&grid.lock[gridCoord]);
}

//
// grid remove
//
bool grid_remove(grid_t & grid, particle_t * p)
{
    int grid_x = grid_coord(p->x);
    int grid_y = grid_coord(p->y);
    int gridCoord = grid_x * grid.size + grid_y;

    // Beginning of critical section
    pthread_mutex_lock(&grid.lock[gridCoord]);

    linkedlist_t * current = grid.grid[gridCoord];

    // No elements?
    if (current == 0)
    {
        // End of critical section
        pthread_mutex_unlock(&grid.lock[gridCoord]);
        return false;
    }

    // Special case for first element
    if (current->value == p)
    {
        grid.grid[gridCoord] = current->next;
        free(current);

        // End of critical section
        pthread_mutex_unlock(&grid.lock[gridCoord]);
        return true;
    }

    linkedlist_t * prev = current;
    current = current->next;

    while(current != 0)
    {
        if (current->value == p)
        {
            prev->next = current->next;
            free(current);

            // End of critical section
            pthread_mutex_unlock(&grid.lock[gridCoord]);
            return true;
        }

        prev = current;
        current = current->next;
    }
    
    // End of critical section
    pthread_mutex_unlock(&grid.lock[gridCoord]);
    return false;
}

//
// clears a grid from values and deallocates any memory from the heap.
//
void grid_clear(grid_t & grid)
{
    for (int i = 0; i < grid.size*grid.size; ++i)
    {
        linkedlist_t * curr = grid.grid[i];
        while(curr != 0)
        {
            linkedlist_t * tmp = curr->next;
            free(curr);
            curr = tmp;
        }
    }
    free(grid.grid);
}

int grid_size(grid_t & grid)
{
    int count = 0;
    for (int i = 0; i < grid.size * grid.size; ++i)
    {
        linkedlist_t * curr = grid.grid[i];
        while(curr != 0)
        {
            count++;
            curr = curr->next;
        }
    }
    return count;
}