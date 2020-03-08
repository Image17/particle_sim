#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;
#include <vector>
//
// particle data structure
//
typedef struct
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  int bx;
  int by;
} particle_t;

typedef struct
{
  std::vector<int> iP;
  std::vector <std::pair<int, int> > blockXY;
} block_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_numblocks( int f );
int get_numblocks();
void set_factor ( int n );
std::pair <int,int> determine_block(double x, double y);
void set_size( int n );
void init_blocks( int n, block_t **blocks, particle_t *p);
void update_blocks( block_t** blocks, particle_t* particles, int n );
void init_particles( int n, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
