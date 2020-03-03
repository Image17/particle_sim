#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <utility>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

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

//
// block data structure
//
typedef struct _block_t
{
  // particle pointer list
  int iP[1000];
  // list of neighboring blocks
  struct _block_t* n_blocks;

  //particle count
  int p_count;
  //index tracker
  int i_track;


  // double x bounds
  double bx_lower;
  double bx_upper;
  // double y bounds
  double by_lower;
  double by_upper;

} block_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p, block_t **blocks );
int get_numblocks();
void update_blocks ( block_t **blocks, particle_t *p, int n );
void init_blocks( int n, block_t **blocks, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p, block_t **blocks, int n  );
std::pair<int, int> determine_block(double x, double y);


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
