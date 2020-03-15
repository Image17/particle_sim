#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "common.h"

double size;
int numblocks;
int factor;

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

double get_size()
{
    return size;
}

void set_factor ( int n )
{
    factor = n;
}

void set_numblocks( int f )
{
    numblocks = ceil(size / (cutoff * f));
}

void update_blocks( block_t** blocks, particle_t* p, int n )
{
  for(int i = 0; i< numblocks; i++)
  {
    for(int j = 0; j < numblocks; j++)
    {
      blocks[i][j].iP.clear();
    }
  }

  for (int i = 0; i < n; i++)
  {
    blocks[p[i].bx][p[i].by].iP.push_back(i);
  }
}

std::pair <int,int> determine_block(double x, double y)
{
  int i = floor(x / (cutoff * factor));
  int j = floor(y / (cutoff * factor));
  return std::make_pair(i, j);
}

std::pair <int,int> determine_thread_block(double x, double y, double block_x_size, double block_y_size)
{
    
   int i = floor(x / block_x_size);
   int j = floor(y / block_y_size);
   //printf("got coords %d , %d \n",i,j);
   return std::make_pair(i, j);
}

void init_blocks( int n, block_t **blocks, particle_t *p)
{
    for (int i = 0; i < numblocks; i++)
    {
        for (int j = 0; j < numblocks; j++)
        {
           blocks[i][j].blockXY.push_back(std::make_pair(i, j));
           if (i-1 >= 0)
           {
              blocks[i][j].blockXY.push_back(std::make_pair(i-1, j));
              if (j-1 >= 0)
              {
                  blocks[i][j].blockXY.push_back(std::make_pair(i-1, j-1));
              }
              if (j+1 < numblocks)
              {
                  blocks[i][j].blockXY.push_back(std::make_pair(i-1, j+1));
              }
           }
           if (i+1 < numblocks)
           {
              blocks[i][j].blockXY.push_back(std::make_pair(i+1, j));
              if (j-1 >= 0)
              {
                  blocks[i][j].blockXY.push_back(std::make_pair(i+1, j-1));
              }
              if (j+1 < numblocks)
              {
                  blocks[i][j].blockXY.push_back(std::make_pair(i+1, j+1));
              }
           }
           if (j-1 >= 0)
           {
              blocks[i][j].blockXY.push_back(std::make_pair(i, j-1));
           }
           if (j+1 < numblocks)
           {
              blocks[i][j].blockXY.push_back(std::make_pair(i, j+1));
           }
        }
    }
}

void clear_out_thread_blocks(std::vector<std::vector<thread_block_t> > thread_blocks, int num_x_blocks, int num_y_blocks)
{
        for (int i = 0; i < num_x_blocks; i++)
    {
        for (int j = 0; j < num_y_blocks; j++)
        {
            thread_blocks[i][j].particles.clear();
            // our sections
            thread_blocks[i][j].top_left_section.clear();
            thread_blocks[i][j].top_section.clear();
            thread_blocks[i][j].top_right_section.clear();
            thread_blocks[i][j].left_section.clear();
            thread_blocks[i][j].right_section.clear();
            thread_blocks[i][j].bottom_left_section.clear();
            thread_blocks[i][j].bottom_section.clear();
            thread_blocks[i][j].bottom_right_section.clear();
            // neighbor sections
            thread_blocks[i][j].top_left_section_neighbor.clear();
            thread_blocks[i][j].top_section_neighbor.clear();
            thread_blocks[i][j].top_right_section_neighbor.clear();
            thread_blocks[i][j].left_section_neighbor.clear();
            thread_blocks[i][j].right_section_neighbor.clear();
            thread_blocks[i][j].bottom_left_section_neighbor.clear();
            thread_blocks[i][j].bottom_section_neighbor.clear();
            thread_blocks[i][j].bottom_right_section_neighbor.clear();
        }
    }
}

void load_particles_into_thread_blocks(int n, std::vector<std::vector<thread_block_t> > thread_blocks, particle_t* particles, double block_x_size, double block_y_size)
{
    // iterate over thread blocks and clear sections
    for (int i =0; i < n; i++)
    {
        particle_t p = particles[i];
        std::pair<int,int> block_coordinates = determine_thread_block(p.x, p.y, block_x_size, block_y_size);
        //printf("checking on stuff.. for %d %d\n",block_coordinates.first,block_coordinates.second);
        thread_block_t ttt = thread_blocks[block_coordinates.first][block_coordinates.second];
        //printf("spill the ttt\n");
        
        // add particle to master list
        thread_blocks[block_coordinates.first][block_coordinates.second].particles.push_back(i);
        
        // determine what border section particle belongs to
        double lower_x_bound = block_coordinates.first * block_x_size;
        double upper_x_bound = (block_coordinates.first + 1) * block_x_size;
        
        double lower_y_bound = block_coordinates.second * block_y_size;
        double upper_y_bound = (block_coordinates.second + 1) * block_y_size;
        
        bool passed_lower_x_bound = (p.x - cutoff) < lower_x_bound;
        bool passed_upper_x_bound = (p.x + cutoff) > upper_x_bound;
                
        bool passed_lower_y_bound = (p.y - cutoff) < lower_y_bound;
        bool passed_upper_y_bound = (p.y + cutoff) > upper_y_bound;
        
        // top right check
        if (passed_upper_x_bound && passed_upper_y_bound)
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].top_right_section.push_back(i);
        } else if (passed_upper_x_bound && passed_lower_y_bound)
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].bottom_right_section.push_back(i);
        } else if (passed_lower_x_bound && passed_upper_y_bound) 
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].top_left_section.push_back(i);
        } else if (passed_lower_x_bound && passed_lower_y_bound) 
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].bottom_left_section.push_back(i);
        }
        
        if (passed_lower_x_bound) 
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].left_section.push_back(i);
        }
        if (passed_upper_x_bound) 
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].right_section.push_back(i);
        }
        if (passed_lower_y_bound)
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].bottom_section.push_back(i);
        }
        if (passed_upper_y_bound)
        {
            thread_blocks[block_coordinates.first][block_coordinates.second].top_section.push_back(i);
        }
    }
}

void init_thread_blocks(int n, std::vector<std::vector<thread_block_t> > thread_blocks, particle_t* particles, int num_x_blocks, int num_y_blocks)
{
    for (int i = 0; i < num_x_blocks; i++)
    {
        for (int j = 0; j < num_y_blocks; j++)
        {
           //TODO we added ourselves as a neighbor, need to account for that
           //blocks[i][j].blockXY.push_back(std::make_pair(i, j));
           
           
           // to our left
           if (i-1 >= 0)
           {
              thread_blocks[i][j].left_section_neighbor = thread_blocks[i-1][j].right_section;
              // bottom left
              if (j-1 >= 0)
              {
                  thread_blocks[i][j].bottom_left_section_neighbor = thread_blocks[i-1][j-1].top_right_section;
              }
              // top left
              if (j+1 < num_y_blocks)
              {
                  thread_blocks[i][j].top_left_section_neighbor = thread_blocks[i-1][j+1].bottom_right_section;
              }
           }
           // to our right
           if (i+1 < num_x_blocks)
           {
              thread_blocks[i][j].right_section_neighbor = thread_blocks[i+1][j].left_section;
              // bottom right
              if (j-1 >= 0)
              {
                  thread_blocks[i][j].bottom_right_section_neighbor = thread_blocks[i+1][j-1].top_left_section;
              }
              // top right
              if (j+1 < num_y_blocks)
              {
                  thread_blocks[i][j].top_right_section_neighbor = thread_blocks[i+1][j+1].bottom_left_section;
              }
           }
           // below 
           if (j-1 >= 0)
           {
              thread_blocks[i][j].bottom_section_neighbor = thread_blocks[i][j-1].top_section;
           }
           // above
           if (j+1 < num_y_blocks)
           {
              thread_blocks[i][j].top_section_neighbor = thread_blocks[i][j+1].bottom_section;
           }
            
            
        }
    }
}

int get_numblocks()
{
  return numblocks;
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p)
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

        std::pair <int,int> blockXY = determine_block(p[i].x, p[i].y);

        p[i].bx = blockXY.first;
        p[i].by = blockXY.second;

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
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
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
void move( particle_t &p )
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
