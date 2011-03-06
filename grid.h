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
	linkedlist_t * * grid;
};

typedef struct grid grid_t;

//
// grid routines
//

int grid_coord(double c);

void grid_init(grid_t & grid, int gridsize);
void grid_add(grid_t & grid, particle_t * particle);
bool grid_remove(grid_t & grid, particle_t * particles);
void grid_clear(grid_t & grid);

#endif