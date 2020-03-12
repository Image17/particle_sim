#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
    std::vector<thread_block_t> flattened_thread_blocks;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    double block_x_size, block_y_size;
    int num_x_blocks, num_y_blocks;
    //thread_block_t** thread_blocks;
    
    numthreads = omp_get_max_threads();
        if (numthreads != 2 && numthreads != 8)
        {
            block_x_size = block_y_size = get_size() / sqrt(numthreads);
            num_x_blocks = num_y_blocks = sqrt(numthreads);
        }
        else
        {
            block_x_size = get_size() / 2;
            block_y_size = get_size() / (numthreads / 2);
            num_x_blocks = 2;
            num_y_blocks = numthreads / 2;
        }
        
        //thread_block_t** thread_blocks = (thread_block_t**)malloc( num_x_blocks * sizeof(thread_block_t) );
        //for (int i = 0; i < num_x_blocks; i++)
        //{
          //thread_blocks[i] = (thread_block_t*) malloc( num_y_blocks * sizeof(thread_block_t) );
        //}
        // load thread_blocks and assign particles into border sections
        //thread_block_t thread_blocks[num_x_blocks][num_y_blocks];
        
        //table = vector<vector<cellValue>> (mLengthX, vector<cellValue>(mLengthY));
        
        
        std::vector<std::vector<thread_block_t>> thread_blocks = std::vector<std::vector<thread_block_t>> (num_x_blocks, std::vector<thread_block_t>(num_y_blocks));
  

        load_particles_into_thread_blocks(n, thread_blocks, particles, block_x_size, block_y_size);
        init_thread_blocks(n, thread_blocks, particles, num_x_blocks, num_y_blocks);
        
        for (int i = 0; i < num_x_blocks; i++)
        {
            for (int j = 0; j < num_y_blocks; j++)
            {
                flattened_thread_blocks.push_back(thread_blocks[i][j]);
            }
        }
        int total_blocks = num_x_blocks * num_y_blocks;

    #pragma omp parallel private(dmin)
    {
        


        numthreads = omp_get_num_threads();
        
        // get number of thread blocks base don numthreads
        // 1:1x1
        // 2:2* split in half
        // 4:2x2
        // 8:8* 2x4
        // 16:4x4
    


    
    // init thread blocks
    for( int step = 0; step < 1000; step++ )
    {
        //printf("===step %d===\n",step);
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for (int i = 0; i < flattened_thread_blocks.size(); i++)
        {
            //printf("for block %d of total blocks %d\n",i,flattened_thread_blocks.size());
            thread_block_t curr_block = flattened_thread_blocks[i];
            // compare every particle in curr_block.particles to others
            for (int j = 0; j < curr_block.particles.size(); j++)
            {
                particles[j].ax = particles[j].ay = 0;
                for (int k = 0; k < curr_block.particles.size(); k++)
                {
                    apply_force( particles[curr_block.particles[j]], particles[curr_block.particles[k]],&dmin,&davg,&navg);
                }
            }
        }
        
		
        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
  
        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	    if (dmin < absmin)
	    { 
	        clear_out_thread_blocks(thread_blocks, num_x_blocks, num_y_blocks);
            load_particles_into_thread_blocks(n, thread_blocks, particles, block_x_size, block_y_size);
            init_thread_blocks(n, thread_blocks, particles, num_x_blocks, num_y_blocks);
	        absmin = dmin; 
	        
	    }
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

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
