/*
________________________________________________
CP2 - PROJEKT: KONEVKTIONS-DIFFUSIONS-GLEICHUNG
P. F. Giesel, M. Neumann
Last Update: 10.07.2019
________________________________________________
Anmerkung : Florian konnte im Verlauf unserer Kooperation nicht mhr bei GitHub pushen!!
*/
#include <bits/stdc++.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <iomanip>
#include <fstream>
#include <chrono>
using namespace std;

bool test_mode  = true;
bool diffusion  = 1;
bool diffusion2 = false; // BOX-MULLER
bool konvektion = 0;
bool discreet   = true;
bool output     = true;

double T1 = 100;
double T2 = 10;

void pt(string dec, double val) // print a value to console (FOR DEBUGGING!!!!)
{
    cout << dec << ": " << setprecision(10)<< val<< endl;
}
double execution_time(auto time) // Timer for debugging
{
    auto end_timer = std::chrono::system_clock::now();
    double diff = chrono::duration_cast<chrono::nanoseconds>(end_timer - time).count();
    diff *= 1e-9;
    cout << "Time taken by program is : "  <<fixed << diff << setprecision(4);
    cout << " sec" << endl;
    return 0;
}
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
    int particle_count  = 4e3;          // number of particles
    int time_step       = 400;          // number of timesteps
    double grid_width   = 20;          // width of the grid (x-koord.)
    double grid_height  = 20;          // height of the grid (y-koord.)
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
    array<array<double,100+1>,100+1> density;
    array<array<double,100+1>,100+1> velocity;
    ofstream out1("position.txt");
    ofstream out2("velocity.txt");
    ofstream out3("density.txt");
    ofstream out4("sum_of_all_particles.txt");
    ofstream out5("set.txt");
    out5<<gridpoints_x << ", "<<gridpoints_y << ", "<<time_step<<", "<<particle_count<<endl;  // n
  // TODO Set distribution of the particles
    int a=grid_height/5; // dis
    int v=3; //velocity for convection
    if(discreet==true)
    {
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
    // SET RANDBEDINGUNG
        // SET TEMP.

    for(int i=1;i<gridpoints_x+1;i++)
    {
        double dens_old_1 = density[1][i];
        density[1][i]=T1;
        particle_count+=T1-dens_old_1;
        for(int m=1;m<T1-dens_old_1+1;m++)
        {
            x_pos.push_back(i);
            y_pos.push_back(1);
            x_vel.push_back(0);
            y_vel.push_back(0);
        }


        double dens_old_2 = density[gridpoints_x+1][i];
        density[gridpoints_x][i]=T2;
        particle_count+=T2-dens_old_2;
        for(int n=1;n<T2-dens_old_2+1;n++)
        {
            x_pos.push_back(i);
            y_pos.push_back(gridpoints_x);
            x_vel.push_back(0);
            y_vel.push_back(0);
        }
    }
        // SET VELOCITY-FIELD
        //Cell-1
    for(int i=a;i<gridpoints_x/2+1;i++)
    {
        velocity[a][i]=v;
        velocity[gridpoints_x-a+1][i]=v;
    }
    for(int j=a;j<gridpoints_y-a+1;j++)
    {
        velocity[j][a]=v;
        velocity[j][gridpoints_y/2]=v;
    }
    //Cell-2
    for(int i=gridpoints_x/2+1;i<gridpoints_x-a+1+1;i++)
    {
        velocity[a][i]=v;
        velocity[gridpoints_x-a+1][i]=v;
    }
    for(int j=a;j<gridpoints_y-a+1;j++)
    {
        velocity[j][gridpoints_y/2+1]=v;
        velocity[j][gridpoints_y-a+1]=v;
    }

    //output t=0
    for (int i=1; i<gridpoints_x+1; i++)
    {
        for (int j=1; j<gridpoints_y+1; j++)
        {
            out3<<i<<", "<<j<<", "<<density[i][j]<<endl;
        }
    }
    for(int g=0;g<gridpoints_x;g++)
    {

    if(v-1>0&&0==1)
    {
        v=v-1;
        for(int i = 0;i<gridpoints_x+1;i++)
        {
            for(int j=0;j<gridpoints_y+1;j++)
            {
                if(velocity[i][j]>v)
                {
                    if(velocity[i+1][j+1]==0)
                    {
                        velocity[i+1][j+1]=v;
                    }
                    if(velocity[i+1][j]==0)
                    {
                        velocity[i+1][j]=v;
                    }
                    if(velocity[i+1][j-1]==0)
                    {
                        velocity[i+1][j-1]=v;
                    }

                    if(velocity[i][j+1]==0)
                    {
                        velocity[i][j+1]=v;
                    }
                    if(velocity[i][j-1]==0)
                    {
                        velocity[i][j-1]=v;
                    }

                    if(velocity[i-1][j+1]==0)
                    {
                        velocity[i-1][j+1]=v;
                    }
                    if(velocity[i-1][j]==0)
                    {
                        velocity[i-1][j]=v;
                    }
                    if(velocity[i-1][j-1]==0)
                    {
                        velocity[i-1][j-1]=v;
                    }
                }
            }
        }
    }
    else{break;}

    }
    //write initila velo-field
    for (int i=1; i<gridpoints_x+1; i++)
    {
        for (int j=1; j<gridpoints_y+1; j++)
        {
            out2<<i<<", "<<j<<", "<<velocity[i][j]<<endl;
        }
    }
  }
    else
    {// NOT DISCRETE
      // SET DISTRIBUT.


      // SET RANDBED.

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
                if(discreet==true)
                {
                    // calculate velocity
                    x_vel[n] = 1; //x_pos_2[n]-x_pos[n];
                    y_vel[n] = 1; //y_pos_2[n]-y_pos[n];

                    //random walk movement with stepsize according to particle velocity
                    //Question: local density could also be taken into account as collision probability rises with higher densities
                    double zufall=dis(gen);
                    if (zufall < 0.25)
                    {
                        if(density[x_pos[n]][y_pos[n]]>0)
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
                    }
                    if (zufall >= 0.25 && zufall < 0.5)
                    {
                        if(density[x_pos[n]][y_pos[n]]>0)
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
                    }
                    if (zufall >= 0.5 && zufall < 0.75)
                    {
                        if(density[x_pos[n]][y_pos[n]]>0)
                        {
                            density[x_pos[n]][y_pos[n]]--;
                            if(y_pos[n]+y_vel[n]==0)
                            {
                                y_pos[n] = y_pos[n]-y_vel[n];
                            }
                            if(y_pos[n]+y_vel[n]==gridpoints_y+1)
                            {
                                y_pos[n] = y_pos[n]-y_vel[n];
                            }
                            else{y_pos[n] = y_pos[n]+y_vel[n];}
                            density[x_pos[n]][y_pos[n]]++;
                        }
                    }
                    if (zufall >= 0.75)
                    {
                        if(density[x_pos[n]][y_pos[n]]>0)
                        {
                            density[x_pos[n]][y_pos[n]]--;
                            if(y_pos[n]-y_vel[n]==0)
                            {
                                y_pos[n] = y_pos[n]+y_vel[n];
                            }
                            if(y_pos[n]-y_vel[n]==gridpoints_y+1)
                            {
                                y_pos[n] = y_pos[n]+y_vel[n];
                            }
                            else{y_pos[n] = y_pos[n]-y_vel[n];}
                            density[x_pos[n]][y_pos[n]]++;
                        }
                    }
                }
                else
                {
                    // NOT DISCRETE

                }
            }


            // convection movement
            if(konvektion==true)
            {
                if(discreet==true)
                { // still just for one grid size - TODO: fit to velocityfield
                    // cell 1
                    for(int i=0;i<2;i++)
                    {
                        if(y_pos[n]>a-1 && y_pos[n]<a+1 && x_pos[n]> a-1 && x_pos[n]< gridpoints_x/2 )
                        {
                            x_pos[n]+=1;
                        }
                        if(y_pos[n]>a && y_pos[n]<gridpoints_y-2 && x_pos[n]> a-1 && x_pos[n]< a-1+2)
                        {
                            y_pos[n]-=1;
                        }
                        if(y_pos[n]<gridpoints_x-a+2 && y_pos[n]>gridpoints_y-a+2-2 && x_pos[n]> a-1+1 && x_pos[n]< gridpoints_x/2 -1+2)
                        {
                            x_pos[n]-=1;
                        }
                        if(y_pos[n]>a-1 && y_pos[n]<gridpoints_y-a+1 && x_pos[n]>gridpoints_x/2 -1 && x_pos[n]<gridpoints_x/2 +1)
                        {
                            y_pos[n]+=1;
                        }
                        //cell-2
                        if(y_pos[n]>a-1 && y_pos[n]<gridpoints_y-a+1 && x_pos[n]> gridpoints_x/2 && x_pos[n]< gridpoints_x/2 +2)
                        {
                            y_pos[n]+=1;
                        }
                        if(y_pos[n]>a-1 && y_pos[n]<a+1 && x_pos[n]<gridpoints_x -a+2 && x_pos[n]> gridpoints_x/2 )
                        {
                            x_pos[n]-=1;
                        }
                        if(y_pos[n]>a && y_pos[n]<gridpoints_y-a+2 && x_pos[n]>gridpoints_x-a && x_pos[n]<gridpoints_x-a+2)
                        {
                            y_pos[n]-=1;
                        }
                        if(y_pos[n]<gridpoints_x-a+2 && y_pos[n]>gridpoints_y-a+2-2 && x_pos[n]<gridpoints_x-a+1 && x_pos[n]> gridpoints_x/2)
                        {
                            x_pos[n]+=1;
                        }
                    }
                }
                else
                {
                    double i_a = 1/a;
                    double i_gx = 1/gridpoints_x;
                    double vgy_a = v*gridpoints_y/a;
                    double i_gy2a = 1/(gridpoints_y*0.5-a);
                    double gra = 1/(gridpoints_y*0.5-a)*gridpoints_y*0.5;
                    if(true==true)
                    {
                        if(x_pos[n]<=gridpoints_x*0.5)
                        {
                            if(y_pos[n]<=a)
                            {
                                x_pos[n] += v*y_pos[n]*i_a;
                            }
                            if(y_pos[n]>a && y_pos[n]<gridpoints_y-a)
                            {
                                x_pos[n] += -v*y_pos[n]*i_gy2a+v*gra;
                                y_pos[n] += v*x_pos[n]*4*i_gx-v;
                            }
                            if(y_pos[n]>=gridpoints_y-a)
                            {
                              x_pos[n] += -v*x_pos[n]*i_a+vgy_a;
                            }
                        }
                        if(x_pos[n]>gridpoints_x*0.5)
                        {
                            if(y_pos[n]<=a)
                            {
                              x_pos[n] += -v*y_pos[n]*i_a;
                            }
                            if(y_pos[n]>a && y_pos[n]<gridpoints_y-a)
                            {
                              x_pos[n] += v*y_pos[n]*i_gy2a-v*gra;
                              y_pos[n] += v*x_pos[n]*4*i_gx-v;
                            }
                            if(y_pos[n]>=gridpoints_y-a)
                            {
                              x_pos[n] += v*x_pos[n]*i_a-vgy_a;
                            }
                        }
                    }
                }

            }
            // write results - pos:
            //if(output==true){
            //out1<<n<<", "<<x_pos[n]<<", "<<y_pos[n]<<endl;}
        }// end of particle loop

if(t==200)
{
    konvektion=true;
    //diffusion=0;
}
//calculate dens new
    if(1==true)
    {

    for(int i=1;i<gridpoints_x+1;i++)
    {
        for(int j=1;j<gridpoints_y+1;j++)
        {
            density[i][j]=0;
        }
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
                    density[j][i]++;
                    break;
                }
            }
        }
    }

    }
       // after all, reset heatcount(particle) to start condition at the top and bottom
           for(int i=1;i<gridpoints_x+1;i++)
            {
                double dens_old_1 = density[1][i];
                density[1][i]=T1;
                if(T1-dens_old_1>0) // set back to T1 by adding
                {
                    //pt("aenderung1",T1-dens_old_1);
                    particle_count+=T1-dens_old_1;
                    for(int m=1;m<T1-dens_old_1+1;m++)
                    {
                        x_pos.push_back(i);
                        y_pos.push_back(1);
                        x_vel.push_back(0);
                        y_vel.push_back(0);
                    }
                }
                else if(T1-dens_old_1<0)
                {
                    particle_count+=T1-dens_old_1;
                    //remove vecs for removed particles
                }


                double dens_old_2 = density[gridpoints_x][i];
                density[gridpoints_x][i]=T2;
                if(T2-dens_old_2>0)
                {
                    //pt("aenderung1",T2-dens_old_2);
                    particle_count+=T2-dens_old_2;
                    for(int n=1;n<T2-dens_old_2+1;n++)
                    {
                        x_pos.push_back(i);
                        y_pos.push_back(gridpoints_x);
                        x_vel.push_back(0);
                        y_vel.push_back(0);
                    }
                }
                else if(T2-dens_old_2<0)
                {
                    particle_count+=T2-dens_old_2;
                    //remove vecs for removed particles

                }
            }
    //TODO write results into .txt files
        if(output==true)
        {
            // velocity:
                //not needed yet. Perhaps in Ph.d.
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
    }//end of timestep loop
    execution_time(timer);
    return 0;
}
//g++ main.cpp -o project_main -fopenmp
