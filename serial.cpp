#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    int blockfactor = read_int( argc, argv, "-f", 1 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    set_factor( blockfactor );
    set_numblocks( blockfactor );
    int numblocks = get_numblocks();
    block_t** blocks = (block_t**)malloc( numblocks * sizeof(block_t));
    for (int i = 0; i < numblocks; i++)
    {
      blocks[i] = (block_t*)malloc( numblocks * sizeof(block_t));
    }
    init_blocks( n, blocks, particles );
    init_particles( n, particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces
        //
        // for( int i = 0; i < n; i++ )
        // {
        //     //particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //     {
        //
        //     }
				//         apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        // }

        for (int i = 0; i < numblocks; i++)
        {
          for (int j = 0; j < numblocks; j++)
          {
            for (int y = 0; y < blocks[i][j].iP.size(); y++)
            {
              int www = blocks[i][j].iP.size();
              //printf("Block: %d %d PCount=%d \n",i,j,www);
              int px = blocks[i][j].iP[y];
              particles[px].ax = particles[px].ay = 0;
              for (int k = 0; k < blocks[i][j].blockXY.size(); k++)
              {
                int bbx = blocks[i][j].blockXY[k].first;
                int bby = blocks[i][j].blockXY[k].second;
                for (int x = 0; x < blocks[bbx][bby].iP.size(); x++)
                {
                  int nx = blocks[bbx][bby].iP[x];
                  apply_force(particles[px], particles[nx], &dmin, &davg, &navg);
                }
              }
            }
          }
        }

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )
            move( particles[i] );

        update_blocks ( blocks, particles, n );

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
