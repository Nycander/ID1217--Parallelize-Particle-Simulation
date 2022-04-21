#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <cstdlib>
#include <pthread.h>
#include <iostream>
using namespace std;

#define NUM_THREADS 4
double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//

void init_particles( int n, particle_t *p ,void *threadarg)
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value)
{

    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;

}


struct thread_data {
   int  thread_id;
   char *message;
   int n;
   particle_t *p ;

};

void *PrintHello(void *threadarg) {
   struct thread_data *my_data;
   my_data = (struct thread_data *) threadarg;
  // sleep(3);

   cout << "Thread ID : " << my_data->thread_id ;
    cout << " Message : " << my_data->message << endl;
   cout << " Ending thread " ;

  pthread_exit(NULL);
}


//void *routine(void *threadarg) {
  // printf(" Test ");
//}

void* init_particles(  void *threadarg )
{
    struct thread_data *args = (struct thread_data *)threadarg;
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)args->n));
    int sy = (args->n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( args->n * sizeof(int) );
    for( int i = 0; i < args->n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < args->n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(args->n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[args->n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        args->p[i].x = size*(1.+(k%sx))/(1+sx);
        args->p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        args->p[i].vx = drand48()*2-1;
        args->p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}
int main (int argc, char* argv[]) {
   pthread_t threads[NUM_THREADS];
   struct thread_data td[NUM_THREADS];
   int rc;
   int i;
   
      cout << "main() start : creating thread[" << 1 << "] next ";
      td[1].thread_id = 1;
      td[1].message = "message";
      td[1].n = 1000;
    //   td[1].p = malloc(sizeof(particle_t)*1000);
      rc = pthread_create(&threads[1], NULL, init_particles, (void *)&td[1]);
      pthread_join(threads[1], NULL);

      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
      cout << "main() start : creating thread[" << 2 << "] next ";
      td[2].thread_id = 2;
      td[2].message = "message";
      td[2].n = 1000;

      rc = pthread_create(&threads[2], NULL, init_particles, NULL );
      pthread_join(threads[2], NULL);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }



      cout << "main() start : creating thread[" << 3 << "] next ";
      td[3].thread_id = 3;
      td[3].message = "This is message";
      td[3].n = 1000;
      rc = pthread_create(&threads[3], NULL, init_particles, (void *)&td[3]);
      pthread_join(threads[3], NULL);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }



      cout << "main() start : creating thread[" << 4 << "] next ";
      td[4].thread_id = 4;
      td[4].message = "message";
      td[4].n = 1000;
      rc = pthread_create(&threads[4], NULL, init_particles, (void *)&td[4]);
      pthread_join(threads[4], NULL);

      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }

      // cout << "creating thread[" << 4 << "] next ";
      // rc = pthread_create(&threads[4], NULL, PrintHello, (void *)4);
      
      // if (rc) {
      //    cout << "Error:unable to create thread," << rc << endl;
      //    exit(-1);
      //    }
   pthread_exit(NULL);
}