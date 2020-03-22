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

typedef struct
{
  std::vector<int> particles;
  int bx_upper;
  int bx_lower;
  int by_upper;
  int by_lower;
  
  std::vector<int> top_left_section;
  std::vector<int> top_left_section_neighbor;
  std::vector<int> top_section;
  std::vector<int> top_section_neighbor;
  std::vector<int> top_right_section;
  std::vector<int> top_right_section_neighbor;
  std::vector<int> left_section;
  std::vector<int> left_section_neighbor;
  std::vector<int> right_section;
  std::vector<int> right_section_neighbor;
  std::vector<int> bottom_left_section;
  std::vector<int> bottom_left_section_neighbor;
  std::vector<int> bottom_section;
  std::vector<int> bottom_section_neighbor;
  std::vector<int> bottom_right_section;
  std::vector<int> bottom_right_section_neighbor;
} thread_block_t;

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
std::vector<std::vector<block_t> > init_blocks_xy( int n, std::vector<std::vector<block_t> > blocks, particle_t *p, int num_x_blocks, int num_y_blocks);
void update_blocks( block_t** blocks, particle_t* particles, int n );
std::vector<std::vector<block_t> > update_blocks_xy( std::vector<std::vector<block_t> > blocks, std::vector<int> p, int n, int num_x_blocks, int num_y_blocks, particle_t* particles );
void init_particles( int n, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );
// openmp
std::pair <int,int> determine_thread_block(double x, double y);
std::vector<std::vector<thread_block_t> > init_thread_blocks(int n, std::vector<std::vector<thread_block_t> > thread_blocks, particle_t* particles, int num_x_blocks, int num_y_blocks);
std::vector<std::vector<thread_block_t> > load_particles_into_thread_blocks(int n, std::vector<std::vector<thread_block_t> > thread_blocks, particle_t* particles, double block_x_size, double block_y_size);
std::vector<std::vector<thread_block_t> > clear_out_thread_blocks(std::vector<std::vector<thread_block_t> > thread_blocks, int num_x_blocks, int num_y_blocks);
double get_size();

void assign_particles_to_blocks(std::vector<int> block_particles, int n, particle_t *p);

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
