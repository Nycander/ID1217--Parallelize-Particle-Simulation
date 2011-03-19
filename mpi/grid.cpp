#include <stdlib.h>
#include <string.h>
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
}

//
// adds a particle pointer to the grid
//
void grid_add(grid_t & grid, particle_t & p)
{
    int gridCoord = grid_coord_flat(grid.size, p.x, p.y);

    linkedlist_t * newElement = (linkedlist_t *) malloc(sizeof(linkedlist));

    newElement->particle_id = p.id;
    newElement->next = grid.grid[gridCoord];

    grid.grid[gridCoord] = newElement;
}

//
// Removes a particle from a grid
//
bool grid_remove(grid_t & grid, particle_t & p, int gridCoord)
{
    if (gridCoord == -1)
        gridCoord = grid_coord_flat(grid.size, p.x, p.y);

    // No elements?
    if (grid.grid[gridCoord] == 0)
    {
        return false;
    }

    linkedlist_t ** nodePointer = &(grid.grid[gridCoord]);
    linkedlist_t * current = grid.grid[gridCoord];

    while(current && (current->particle_id != p.id))
    {
        nodePointer = &(current->next);
        current = current->next;
    }

    if (current)
    {
        *nodePointer = current->next;
        free(current);
    }

    return !!current;
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

void grid_clear_row(grid_t & grid, int r)
{
    for (int c = 0; c < grid.size; ++c)
    {
        int i = r * grid.size + c;
        linkedlist_t * curr = grid.grid[i];
        while(curr != 0)
        {
            linkedlist_t * tmp = curr->next;
            free(curr);
            curr = tmp;
        }
        memset(grid.grid + i, 0, sizeof(linkedlist_t*));
    }
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
/*
int grid_print(grid_t & grid)
{
    for (int i = 0; i < grid.size * grid.size; ++i)
    {
        linkedlist_t * curr = grid.grid[i];
        while(curr != 0)
        {
            printf("x:%d y:%d\n", curr->particle)        
            curr = curr->next;
        }
    }
    return count;
}
*/