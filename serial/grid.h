#ifndef __GRID_H_
#define __GRID_H_

#include "common.h"

struct linkedlist
{
	linkedlist * next;
	particle_t * value;
};

typedef struct linkedlist linkedlist_t;

struct grid
{
	int size;
	linkedlist_t ** grid;
};

typedef struct grid grid_t;

//
// grid routines
//

void grid_init(grid_t & grid, int gridsize);
void grid_add(grid_t & grid, particle_t * particle);
bool grid_remove(grid_t & grid, particle_t * p, int gridCoord = -1);
void grid_clear(grid_t & grid);
int  grid_size(grid_t & grid);


//
// Calculate the grid coordinate from a real coordinate
//
inline static int grid_coord(double c)
{
    return (int)floor(c / cutoff);
}
inline static int grid_coord_flat(int size, double x, double y)
{
    return grid_coord(x) * size + grid_coord(y);
}

#endif