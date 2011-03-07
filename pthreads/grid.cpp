#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include "grid.h"


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
    int gridCoord = grid_coord_flat(grid.size, p->x, p->y);

    linkedlist_t * newElement = (linkedlist_t *) malloc(sizeof(linkedlist));
    newElement->value = p;

    // Beginning of critical section
    pthread_mutex_lock(&grid.lock[gridCoord]);
    newElement->next = grid.grid[gridCoord];

    grid.grid[gridCoord] = newElement;
    // End of critical section
    pthread_mutex_unlock(&grid.lock[gridCoord]);
}

//
// Removes a particle from a grid
//
bool grid_remove(grid_t & grid, particle_t * p)
{
    int gridCoord = grid_coord_flat(grid.size, p->x, p->y);

    // No elements?
    if (grid.grid[gridCoord] == 0)
    {
        return false;
    }

    // Special case for first element
    pthread_mutex_lock(&grid.lock[gridCoord]); // Beginning of critical section
    if (grid.grid[gridCoord]->value == p)
    {
        linkedlist_t * tmp = grid.grid[gridCoord];
        grid.grid[gridCoord] = tmp->next;
        free(tmp);
        
        pthread_mutex_unlock(&grid.lock[gridCoord]); // End of critical section
        return true;
    }

    linkedlist_t * prev = grid.grid[gridCoord];
    linkedlist_t * current = prev->next;
    while(current != 0)
    {
        if (current->value == p)
        {
            prev->next = current->next;
            free(current);

            pthread_mutex_unlock(&grid.lock[gridCoord]); // End of critical section
            return true;
        }

        prev = current;
        current = current->next;
    }
    pthread_mutex_unlock(&grid.lock[gridCoord]); // End of critical section
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