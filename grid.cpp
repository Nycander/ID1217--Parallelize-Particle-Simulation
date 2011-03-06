#include <stdlib.h>
#include <string.h>

#include "grid.h"


//
// Calculate the grid coordinate from a real coordinate
//
int grid_coord(double c)
{
	return (int)(c / cutoff);
}


//
// initialize grid and fill it with particles
// 
void grid_init(grid_t & grid, int size)
{
	grid.size = size;
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
void grid_add(grid_t & grid, particle_t * p)
{
	int gridx = grid_coord(p->x);
	int gridy = grid_coord(p->y);

	int gridCoord = gridx * grid.size + gridy;
	
	linkedlist_t * tmp = grid.grid[gridCoord];

	linkedlist_t * newElement = (linkedlist_t *) malloc(sizeof(linkedlist));
	newElement->next = tmp;
	newElement->value = p;

	grid.grid[gridCoord] = newElement;
}

//
// grid move
//
bool grid_remove(grid_t & grid, particle_t * p)
{
	int grid_x = grid_coord(p->x);
	int grid_y = grid_coord(p->y);
	int gridCoord = grid_x * grid.size + grid_y;

	// Find p in linkedlist
	linkedlist_t * current = grid.grid[gridCoord];

	// No elements?
	if (current == 0)
	{
		return false;
	}
	
	// Special case for size = 1
	if (current->next == 0)
	{
		if (current->value != p)
		{
			return false;
		}

		grid.grid[gridCoord] = 0;
		free(current);
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
			return true;
		}

		prev = current;
		current = current->next;
	}
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