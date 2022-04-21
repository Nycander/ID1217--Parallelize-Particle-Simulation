#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include <pthread.h>
using namespace std;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

// calculate particle's bin number
int binNum(particle_t &p, int bpr) 
{
    return ( floor(p.x/cutoff) + bpr*floor(p.y/cutoff) );
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    // create spatial bins (of size cutoff by cutoff)
    double size = sqrt( density*n );
    int bpr = ceil(size/cutoff);
    int numbins = bpr*bpr;
    vector<particle_t*> *bins = new vector<particle_t*>[numbins];

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {

      // clear bins at each time step
      for (int m = 0; m < numbins; m++)
	bins[m].clear();
    
      // place particles in bins
      for (int i = 0; i < n; i++) 
	bins[binNum(particles[i],bpr)].push_back(particles + i);

      //
      //  compute forces
      //
      for( int p = 0; p < n; p++ )
        {
          particles[p].ax = particles[p].ay = 0;
          
	  // find current particle's bin, handle boundaries
	  int cbin = binNum( particles[p], bpr );
	  int lowi = -1, highi = 1, lowj = -1, highj = 1;
	  if (cbin < bpr)
	    lowj = 0;
	  if (cbin % bpr == 0)
	    lowi = 0;
	  if (cbin % bpr == (bpr-1))
	    highi = 0;
	  if (cbin >= bpr*(bpr-1))
	    highj = 0;

	  // apply nearby forces
	  for (int i = lowi; i <= highi; i++)
	    for (int j = lowj; j <= highj; j++)   
	      {
		int nbin = cbin + i + bpr*j;
		for (int k = 0; k < bins[nbin].size(); k++ )
		  apply_force( particles[p], *bins[nbin][k] );
	      }
        }
        
        //
        //  move particles
        //
        for( int p = 0; p < n; p++ ) 
            move( particles[p] );
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    delete [] bins;
    if( fsave )
        fclose( fsave );
    
    return 0;
}
