#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <utility>
#include <iostream>

using namespace std;

double size;
int numblocks;

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
    numblocks = (int)ceil(size/cutoff);
}

int get_numblocks()
{
  return numblocks;
}

void init_blocks( int n, block_t **blocks, particle_t *p)
{

    for (int i = 0; i < numblocks; i++)
    {
        for (int j = 0; j < numblocks; j++)
        {

           //set bounds
           // set upper and lower bounds for x coordinates
           blocks[i][j].bx_lower = cutoff * (i);
           blocks[i][j].bx_upper = cutoff * (i + 1);

           // set upper and lower bounds for y coordinates
           blocks[i][j].by_lower = cutoff * (j);
           blocks[i][j].by_upper = cutoff * (j + 1);


           if (i == 0 || i == (numblocks-1) || j == 0 || j == numblocks-1)
           {
              if ((i == 0 && j == 0) || (i == 0 && j == (numblocks-1)) || (i == (numblocks-1) && j == 0) || (i == (numblocks-1) && j == (numblocks-1)))
              {
                  //Corners
                  blocks[i][j].n_blocks = (block_t *)malloc(3 * sizeof(block_t));
              }
              else
              {
                  //Sides
                  blocks[i][j].n_blocks = (block_t *)malloc(5 * sizeof(block_t));
              }
           }
           else
           {
              //Inner
              blocks[i][j].n_blocks = (block_t *)malloc(8 * sizeof(block_t));
           }
           //set neighbor blocks
           int k = 0;
           if (i-1 >= 0)
           {
              blocks[i][j].n_blocks[k] = blocks[i-1][j];
              k++;
              if (j-1 >= 0)
              {
                  blocks[i][j].n_blocks[k] = blocks[i-1][j-1];
                  k++;
              }
              if (j+1 < numblocks)
              {
                  blocks[i][j].n_blocks[k] = blocks[i-1][j+1];
                  k++;
              }
           }
           if (i+1 < numblocks)
           {
              blocks[i][j].n_blocks[k] = blocks[i+1][j];
              k++;
              if (j-1 >= 0)
              {
                  blocks[i][j].n_blocks[k] = blocks[i+1][j-1];
                  k++;
              }
              if (j+1 < numblocks)
              {
                  blocks[i][j].n_blocks[k] = blocks[i+1][j+1];
                  k++;
              }
           }
           if (j-1 >= 0)
           {
              blocks[i][j].n_blocks[k] = blocks[i][j-1];
              k++;
           }
           if (j+1 < numblocks)
           {
              blocks[i][j].n_blocks[k] = blocks[i][j+1];
              k++;
           }

           //set inital particle lists
           //if within bounds, add to list
           blocks[i][j].p_count = 0;
           //load_block(blocks[i][j], p, n);
        }
    }
}

void update_blocks ( block_t **blocks, particle_t *p, int n )
{
    //after move and init update block particle lists
    //reset all block p lists with new arrays

    for (int i = 0; i < numblocks; i++)
    {
        for (int j = 0; j < numblocks; j++)
        {
            blocks[i][j].i_track = 0;
            //blocks[i][j].iP = (int *) malloc (blocks[i][j].p_count * sizeof(int));
            //blocks[i][j].iP = new int[n];
        }
    }

    //Loop through particles and assign them to their n_blocks

    for (int i = 0; i < n; i++)
    {
        //cout << "adding particle at " << i << "to our x y " << p[i].bx << p[i].by << endl;
        int k = blocks[p[i].bx][p[i].by].i_track;
        blocks[p[i].bx][p[i].by].iP[k] = i;
        blocks[p[i].bx][p[i].by].i_track++;
    }
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p, block_t **blocks )
{
    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));

    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
    //printf("for i = %d\n",i);
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
    
        std::pair<int, int> blockXY = determine_block(p[i].x, p[i].y);
    //printf("for our coords %d, %d we are going to assign to block %d, %d\n",p[i].x,p[i].y,blockXY.first,blockXY.second);
    //printf("wtf %d\n", blockXY.first);
        blocks[blockXY.first][blockXY.second].p_count++;
    //cout << "pair x: " << blockXY.first << endl;
    //cout << "pair y: " << blockXY.second << endl;
    //cout << "our p count " << blocks[blockXY.first][blockXY.second].p_count << endl;


        p[i].bx = blockXY.first;
        p[i].by = blockXY.second;
    //p[i].bx = p[i].x;
    //p[i].by = p[i].y;
        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
        //printf(" %d ", i);
    }

    update_blocks(blocks, p, n);

    free( shuffle );
}


std::pair<int, int> determine_block(double x, double y)
{
  int i = (int)floor(x / cutoff);
  int j = (int)floor(y  / cutoff);
  std::pair<int, int> blockXY = std::make_pair(i, j);
  //printf("determine coords %d, %d we are translating into %d, %d and then into pair %d, %d\n",x,y,i,j,blockXY.first,blockXY.second);
  return blockXY;
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    //cout << "in apply force...." << endl;
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
  //printf("calculations: %d * %d = %d\n",dx,dy,r2);
    if( r2 > cutoff*cutoff ) {
        //printf("too far, skipping");
    return;
  }
  if (r2 != 0)
        {
     if (r2/(cutoff*cutoff) < *dmin * (*dmin))
        //cout << "changing dmin" << endl;
        *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }

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
void move( particle_t &p, block_t **blocks, int n )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    int oldBX = p.bx;
    int oldBY = p.by;

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
    std::pair<int, int> blockXY = determine_block(p.x, p.y);
    if (oldBX != blockXY.first || oldBY != blockXY.second)
    {
        p.by = blockXY.second;
        p.bx = blockXY.first;
        if (blocks[oldBX][oldBY].p_count > 0)
          blocks[oldBX][oldBY].p_count--;
        if (blocks[blockXY.first][blockXY.second].p_count < n)
          blocks[blockXY.first][blockXY.second].p_count++;
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

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
