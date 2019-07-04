//Du kannst auch ganz von Vorne anfangen, wenn dir das besser liegt Marvin
//Hier ist bis jetzt nur der unfertige Zufallsspaziergang drin
#include <fstream>
#include <ofstream>
#include <vector>
#include <random>
using namespace std;

int main()
{
//integration of random number generator
  mt19937 gen;
  gen.seed(45);

  uniform_real_distribution<double> dis(0,1);

//predefining variables and vecors
  const int quantity=1000; //number of particles
  const int steps=100; //number of timesteps
  vector<int> x_position(quantity,0);
  vector<int> y_position(quantity,0);
  vector<int> x1_position(quantity,0);
  vector<int> y1position(quantity,0);

  vector<int> speed(quantity,0);

//loop over all steps
  for (int t=0; t<steps; t++)
  {
//loop over all particles moving particles for each step
    for (int n=0; n<quantity; n++)
    {
//TODO calculating current particle speed. particle density and stepsize for random walk
//Problem: calculating particle speed by v=sqrt((dx/dt)²+(dy/dt)²) would produce non-integers
      speed[n]=;
//writing current particle positions into new vector
      x1_position[n] = x_position[n];
      y1_position[n] = y_position[n];
//random walk movement with stepsize according to particle speed
//Question: local density could also be taken into account as collision probability rises with higher densities
      double zufall=dis(gen);
      if (zufall < 0.25){
	x_position[n] = x_position[n]+speed[n];}
      if (zufall >= 0.25 && zufall < 0.5){
	x_position[n] = x_position[n]-speed[n];}
      if (zufall >= 0.5 && zufall < 0.75){
	y_position[n] = y_position[n]+speed[n];}
      if (zufall >= 0.75){
	y_position[n] = y_position[n]-speed[n];}
//TODO convection movement
    }
  }
//write results into txt-files
  ofstream out("ort.txt");
  for (int i=0; i<quantity; i++)
  {
    out<<x_position[i]<<"\t"<<y_position[i]<<endl;
  }
}
