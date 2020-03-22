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
    printf("in main\n");
    int navg,nabsavg=0,numthreads,numblocks;
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


        std::vector<std::vector<thread_block_t> > thread_blocks = std::vector<std::vector<thread_block_t> > (num_x_blocks, std::vector<thread_block_t>(num_y_blocks));


        thread_blocks = load_particles_into_thread_blocks(n, thread_blocks, particles, block_x_size, block_y_size);
        thread_blocks = init_thread_blocks(n, thread_blocks, particles, num_x_blocks, num_y_blocks);

        for (int i = 0; i < num_x_blocks; i++)
        {
            for (int j = 0; j < num_y_blocks; j++)
            {
                flattened_thread_blocks.push_back(thread_blocks[i][j]);
            }
        }
        int total_blocks = num_x_blocks * num_y_blocks;
        
    printf("before parallel\n");
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
    for( int step = 0; step < 1; step++ )
    {
        printf("===step %d===\n",step);
        navg = 0;
        davg = 0.0;
	    dmin = 1.0;
        //
        //  compute all forces
        //
            thread_block_t curr_block = flattened_thread_blocks[omp_get_thread_num()];
            // compare every particle in curr_block.particles to others
            //printf("for step %d and thread %d we have %d particles \n", step, omp_get_thread_num(), curr_block.particles.size());
            set_factor( 1 );
            set_numblocks( 1 );
            int numblocks = get_numblocks();
            int xblocks = ceil(block_x_size / .01);
            int yblocks = ceil(block_y_size / .01);

            
             std::vector<std::vector<block_t> > blocks = std::vector<std::vector<block_t> > (xblocks, std::vector<block_t>(yblocks));
            
            printf("finished instantiating blocks\n");
            
            std::vector<particle_t> my_p = std::vector<block_t>(curr_block.particles.size());
            int actualpsize=0;
            for (int i = 0; i < curr_block.particles.size(); i++)
            {
                my_p[i] = particles[curr_block.particles[i]];
                actualpsize++;
            }
             printf("finished instantiating particles, we have %d particles\n",actualpsize);
             
             printf("starting init_blocks\n");
            init_blocks_xy( curr_block.particles.size(), blocks, particles, xblocks, yblocks );
            printf("finished init_blocks\n");
            
            printf("starting update\n");
            update_blocks_xy( blocks, my_p, actualpsize, xblocks, yblocks );
            printf("ending update\n");

            
            printf("beginning on my blocks..\n");
            for (int i = 0; i < xblocks; i++)
            {
              for (int j = 0; j < yblocks; j++)
              {
                  if (blocks[i][j].iP.size() > 0){
                    printf("the current block for %d %d has %d particles\n",i,j,blocks[i][j].iP.size());
                  }
                for (int y = 0; y < blocks[i][j].iP.size(); y++)
                {
                  int px = blocks[i][j].iP[y];
                  particles[px].ax = particles[px].ay = 0;
                  for (int k = 0; k < blocks[i][j].blockXY.size(); k++)
                  {
                    int bbx = blocks[i][j].blockXY[k].first;
                    int bby = blocks[i][j].blockXY[k].second;
                    for (int x = 0; x < blocks[bbx][bby].iP.size(); x++)
                    {
                        printf("applying force\"");
                      int nx = blocks[bbx][bby].iP[x];
                      apply_force(particles[px], particles[nx], &dmin, &davg, &navg);
                    }
                  }
                }
              }
            }
            printf("out of apply force loop\n");
            


          //  printf("for step %d and thread %d we had %d interactions \n", step, omp_get_thread_num(), interactios);
         //#pragma omp for reduction (+:navg) reduction(+:davg)
         printf("starting right\n");
         for (int ri = 0; ri < curr_block.right_section.size(); ri++)
            {
                for (int rin = 0; rin < curr_block.right_section_neighbor.size(); rin++)
                {
                    apply_force( particles[curr_block.right_section[ri]], particles[curr_block.right_section_neighbor[rin]],&dmin,&davg,&navg);
                }
                for (int rtin = 0; rtin < curr_block.top_right_section_neighbor.size(); rtin++)
                {
                    apply_force( particles[curr_block.right_section[ri]], particles[curr_block.top_right_section_neighbor[rtin]],&dmin,&davg,&navg);
                }
                for (int rbin = 0; rbin < curr_block.bottom_right_section_neighbor.size(); rbin++)
                {
                    apply_force( particles[curr_block.right_section[ri]], particles[curr_block.bottom_right_section_neighbor[rbin]],&dmin,&davg,&navg);
                }
            }
         //#pragma omp for reduction (+:navg) reduction(+:davg)
         printf("starting left\n");
         for (int li = 0; li < curr_block.left_section.size(); li++)
            {
                for (int lin = 0; lin < curr_block.left_section_neighbor.size(); lin++)
                {
                    apply_force( particles[curr_block.left_section[li]], particles[curr_block.left_section_neighbor[lin]],&dmin,&davg,&navg);
                }
                for (int ltin = 0; ltin < curr_block.top_left_section_neighbor.size(); ltin++)
                {
                    apply_force( particles[curr_block.left_section[li]], particles[curr_block.top_left_section_neighbor[ltin]],&dmin,&davg,&navg);
                }
                for (int lbin = 0; lbin < curr_block.bottom_left_section_neighbor.size(); lbin++)
                {
                    apply_force( particles[curr_block.left_section[li]], particles[curr_block.bottom_left_section_neighbor[lbin]],&dmin,&davg,&navg);
                }
            }
         //#pragma omp for reduction (+:navg) reduction(+:davg)
         printf("starting top\n");
         for (int li = 0; li < curr_block.top_section.size(); li++)
            {
                for (int lin = 0; lin < curr_block.top_section_neighbor.size(); lin++)
                {
                    apply_force( particles[curr_block.top_section[li]], particles[curr_block.top_section_neighbor[lin]],&dmin,&davg,&navg);
                }
                for (int ltin = 0; ltin < curr_block.top_right_section_neighbor.size(); ltin++)
                {
                    apply_force( particles[curr_block.top_section[li]], particles[curr_block.top_right_section_neighbor[ltin]],&dmin,&davg,&navg);
                }
                for (int lbin = 0; lbin < curr_block.top_left_section_neighbor.size(); lbin++)
                {
                    apply_force( particles[curr_block.top_section[li]], particles[curr_block.top_left_section_neighbor[lbin]],&dmin,&davg,&navg);
                }
            }
         //#pragma omp for reduction (+:navg) reduction(+:davg)
         printf("starting bottom\n");
         for (int li = 0; li < curr_block.bottom_section.size(); li++)
            {
                for (int lin = 0; lin < curr_block.bottom_section_neighbor.size(); lin++)
                {
                    apply_force( particles[curr_block.bottom_section[li]], particles[curr_block.bottom_section_neighbor[lin]],&dmin,&davg,&navg);
                }
                for (int ltin = 0; ltin < curr_block.bottom_right_section_neighbor.size(); ltin++)
                {
                    apply_force( particles[curr_block.bottom_section[li]], particles[curr_block.bottom_right_section_neighbor[ltin]],&dmin,&davg,&navg);
                }
                for (int lbin = 0; lbin < curr_block.bottom_left_section_neighbor.size(); lbin++)
                {
                    apply_force( particles[curr_block.bottom_section[li]], particles[curr_block.bottom_left_section_neighbor[lbin]],&dmin,&davg,&navg);
                }
            }
        printf("out of neighbor loop\n");
        #pragma omp barrier


        //
        //  move particles
        //
        #pragma omp for
        for( int i = 0; i < n; i++ )
            move( particles[i] );
        printf("done move\n");

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

          #pragma omp single
          {
            thread_blocks = clear_out_thread_blocks(thread_blocks, num_x_blocks, num_y_blocks);
            thread_blocks = load_particles_into_thread_blocks(n, thread_blocks, particles, block_x_size, block_y_size);
            thread_blocks = init_thread_blocks(n, thread_blocks, particles, num_x_blocks, num_y_blocks);
            flattened_thread_blocks.clear();
           for (int i = 0; i < num_x_blocks; i++)
            {
                for (int j = 0; j < num_y_blocks; j++)
                {
                    flattened_thread_blocks.push_back(thread_blocks[i][j]);
                }
            }
          }
          #pragma omp critical
          {
    	    if (dmin < absmin)
    	    {

    	        absmin = dmin;

    	    }
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
