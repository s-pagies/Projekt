/*
________________________________________________
CP2 - PROJEKT: KONEVKTIONS-DIFFUSIONS-GLEICHUNG
P. F. Giesel, M. Neumann
Last Update: 06.07.2019
________________________________________________
ANNAHMEN:
- Homogene Masse (Besteht nur aus einem Stoff, anfangs überall gleich bzgl. Dichte<->Temperatur)
- 2 Dimensional
- Umgebungsdruck konstant       = 1. GP
- Gravitationkonstante konstant = 9.81m/s^2
-
-
*/
#include <bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <fstream>
#include <chrono>
using namespace std;

bool test_mode = true;
bool diffusion = true;
bool konvektion= true;
bool output    = true;
void pt(string dec, double val) // print a value to console (FOR DEBUGGING!!!!)
{
    cout << dec << ": " << val << endl;
}
// Timer for debugging
//auto timer = chrono::system_clock::now();
double execution_time(auto time)
{
    auto end_timer = std::chrono::system_clock::now();
    double diff = chrono::duration_cast<chrono::nanoseconds>(end_timer - time).count();
    diff *= 1e-9;
    cout << "Time taken by program is : "  <<fixed << diff << setprecision(4);
    cout << " sec" << endl;
    return 0;
}
double execution_time2(auto time)
{
    auto end_timer = std::chrono::system_clock::now();
    chrono::duration<double> diff = end_timer-time;
    time_t end_time = chrono::system_clock::to_time_t(end_timer);
    pt("exeution_time in (s)", diff.count());
    return 0;
}

float g = 9.81; // gravity-const.
double m = 1.0;
double F_g = g*m;
//#####################################################################
int main()
{
    auto timer = chrono::system_clock::now();
    // integration of random number generator
    double seed = 42.0;
    mt19937 gen;
    gen.seed(seed);
    uniform_real_distribution<double> dis(0,1);

    //predefining variables, vectors for particles and grid
    int particle_count  = 1e5;          // number of particles
    int time_step       = 10;          // number of timesteps
    double grid_width   = 100;          // width of the grid (x-koord.)
    double grid_height  = 100;          // height of the grid (y-koord.)
    double grid_x_step  = 1;            // stepwidth for x
    double grid_y_step  = 1;            // stepwidth for y
    double gridpoints_x = grid_width/grid_x_step;
    double gridpoints_y = grid_height/grid_y_step;
    double gridpoint_count = gridpoints_x*gridpoints_y;
    double particles_per_gridpoint = particle_count/gridpoint_count;

    if(particles_per_gridpoint<1)
    {
        pt("error: particles_per_gridpoint is < 1 ",0);
        return 0;
    }

    vector<double> x_pos(particle_count+1,0.0);
    vector<double> y_pos(particle_count+1,0.0);
    vector<double> x_pos_2(particle_count+1,0.0);
    vector<double> y_pos_2(particle_count+1,0.0);
    vector<double> x_vel(particle_count+1,0.0);
    vector<double> y_vel(particle_count+1,0.0);
    vector<double> velocity(particle_count+1,0.0);
    array<array<double,100+1>,100+1> density;
    ofstream out1("position.txt");
    ofstream out2("velocity.txt");
    ofstream out3("density.txt");
    ofstream out4("sum_of_all_particles.txt");
  // TODO Set distribution of the particles

  // Equal distribution:
    int abstand = 0;
    for(int i=1;i<gridpoints_x+1;i++)
    {
        for(int j=1;j<gridpoints_y+1;j++)
        {
            // set all density to 0
            density[i][j]=0;
            // set position of the particles
            for(int k=1;k<particles_per_gridpoint+1;k++)
            {
                x_pos[k+particles_per_gridpoint*i+particles_per_gridpoint*j-2*particles_per_gridpoint+(gridpoints_y-1)*abstand]=i;
                y_pos[k+particles_per_gridpoint*i+particles_per_gridpoint*j-2*particles_per_gridpoint+(gridpoints_y-1)*abstand]=j;
            }
        }
        abstand+=particles_per_gridpoint;
    }

    // get density from position
    for(int n =1;n<particle_count+1;n++)
    {
        for(int i=1;i<gridpoints_x+1;i++)
        {
            for(int j=1;j<gridpoints_y+1;j++)
            {
                if(x_pos[n]==i && y_pos[n]==j)
                {
                    density[i][j]++;
                    break;
                }
            }
        }
    }

    if(test_mode==true)
    {
        pt("particle_count", particle_count);
        pt("time_steps", time_step);
        pt("gridpoints_x", gridpoints_x);
        pt("gridpoints_y", gridpoints_y);
        pt("particles_per_gridpoint", particles_per_gridpoint);
        cout << "######################################################" << '\n'<<endl;
    }

    //loop over all time steps
    for (int t=1; t<time_step+1; t++)
    {
        pt("t", t);
        execution_time(timer);
        //loop over all particles moving particles for each step
        for (int n=1; n<particle_count+1; n++)
        {
            if(test_mode==true)
            {
                if(n==1)
                {
                    pt("teilchen", n);
                }
                if(n==particle_count/10)
                {
                    pt("teilchen", n);
                }
                if(n==particle_count*25/100)
                {
                    pt("teilchen", n);
                }
                if(n==particle_count*5/10)
                {
                    pt("teilchen", n);
                }
                if(n==particle_count*75/100)
                {
                    pt("teilchen", n);
                }
                if(n==particle_count)
                {
                    pt("teilchen", n);
                }
            }
            //Diffusion
            if(diffusion==true)
            {
                //TODO calculating current particle velocity, particle density and stepsize for random walk
                //Problem: calculating particle velocity by v=sqrt((dx/dt)²+(dy/dt)²) would produce non-integers
                velocity[n]=1;
                //writing current particle positions into new vector to calculate velocity in next timestep
                x_pos_2[n] = x_pos[n];
                y_pos_2[n] = y_pos[n];

                // calculate velocity
                x_vel[n] = 1; //x_pos_2[n]-x_pos[n];
                y_vel[n] = 1; //y_pos_2[n]-y_pos[n];

                //random walk movement with stepsize according to particle velocity
                //Question: local density could also be taken into account as collision probability rises with higher densities
                double zufall=dis(gen);
                if (zufall < 0.25)
                    {
                        density[x_pos[n]][y_pos[n]]--;
                        if(x_pos[n]+x_vel[n]==0)
                        {
                            x_pos[n] = x_pos[n]-x_vel[n];
                        }
                        if(x_pos[n]+x_vel[n]==gridpoints_x+1)
                        {
                            x_pos[n] = x_pos[n]-x_vel[n];
                        }
                        else{x_pos[n] = x_pos[n]+x_vel[n];}
                        density[x_pos[n]][y_pos[n]]++;
                    }
                if (zufall >= 0.25 && zufall < 0.5)
                    {
                        density[x_pos[n]][y_pos[n]]--;
                        if(x_pos[n]-x_vel[n]==0)
                        {
                            x_pos[n] = x_pos[n]+x_vel[n];
                        }
                        if(x_pos[n]-x_vel[n]==gridpoints_x+1)
                        {
                            x_pos[n] = x_pos[n]+x_vel[n];
                        }
                        else{x_pos[n] = x_pos[n]-x_vel[n];}
                        density[x_pos[n]][y_pos[n]]++;
                    }
                if (zufall >= 0.5 && zufall < 0.75)
                {
                        density[x_pos[n]][y_pos[n]]--;
                        if(y_pos[n]+y_vel[n]==0)
                        {
                            y_pos[n] = y_pos[n]-y_vel[n];
                        }
                        if(y_pos[n]+y_vel[n]==gridpoints_y+1)
                        {
                            x_pos[n] = y_pos[n]-y_vel[n];
                        }
                        else{y_pos[n] = y_pos[n]+y_vel[n];}
                        density[x_pos[n]][y_pos[n]]++;
                }
                if (zufall >= 0.75)
                {
                        density[x_pos[n]][y_pos[n]]--;
                        if(y_pos[n]-y_vel[n]==0)
                        {
                            y_pos[n] = y_pos[n]+y_vel[n];
                        }
                        if(y_pos[n]-y_vel[n]==gridpoints_y+1)
                        {
                            x_pos[n] = y_pos[n]+y_vel[n];
                        }
                        else{y_pos[n] = y_pos[n]-y_vel[n];}
                        density[x_pos[n]][y_pos[n]]++;
                }
            }
                        //TODO convection movement
            if(konvektion==true)// konvecktionSTERM = -v \nabla(c)
            {
                double c        = density[x_pos[n]][y_pos[n]];
                if(x_pos[n]==1)// if particle is at the bottom of the room, don't check the density of the point below
                {
                    double dc_x_p   = density[x_pos[n]+1][y_pos[n]]-c/grid_x_step;
                    //double dc_x_n   = NULL;
                    double dc_x     = dc_x_p;
                }
                if(x_pos[n]==grid_width) // if particle is on top of the room, don't check the density of the point above
                {
                    //double dc_x_p   = NULL;
                    double dc_x_n   = density[x_pos[n]-1][y_pos[n]]-c/grid_x_step;
                    double dc_x     = dc_x_n;
                }
                else
                {
                    double dc_x_p   = density[x_pos[n]+1][y_pos[n]]-c/grid_x_step;
                    double dc_x_n   = density[x_pos[n]-1][y_pos[n]]-c/grid_x_step;
                    double dc_x     = 1/2*(dc_x_p-dc_x_n);
                }
                double dc_y_p   = density[x_pos[n]][y_pos[n]+1]-c/grid_y_step;
                double dc_y_n   = density[x_pos[n]-1][y_pos[n]-1]-c/grid_y_step;
                double dc_y     = 1/2*(dc_y_p-dc_y_n);
                //dc[n]
            }
            // write results - pos:
            if(output==true){
            out1<<n<<", "<<x_pos[n]<<", "<<y_pos[n]<<endl;}
        }


    //TODO write results into .txt files
        if(output==true)
        {
            // velocity:

            // density:
            int sum=0;
            for (int i=1; i<gridpoints_x+1; i++)
            {
                for (int j=1; j<gridpoints_y+1; j++)
                {
                    if(test_mode==true)
                    {
                        sum+=density[i][j];
                    }
                    out3<<i<<", "<<j<<", "<<density[i][j]<<endl;
                }
            }
            out4<<t<<", "<< sum<<endl;
        }

    }
    execution_time(timer);
    return 0;
}
