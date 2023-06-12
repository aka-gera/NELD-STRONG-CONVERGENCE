#include "input.hpp"
#include "matrix.hpp"

class particle
{
public:
  //---------------------------------
  //  All the masses are set to 1 !
  //---------------------------------

  Matrix link;

  Matrix q;    // positions
  Matrix p;    // momenta
  Matrix prel; // relative momentum
  Matrix RF;   // noise term
  Matrix RF2;  // noise term
  Matrix RF3;  // noise term
  Matrix RF4;  // noise term
  Matrix f;    // forces on all particles

  Matrix sigma, pressure;

  double energy;           // total energy
  double potential_energy; // potential energy
  double total_momentum_x;
  double total_momentum_y;
  double total_momentum_z;
  double background_error;

  double mindist;

  double nb_birth, nb_death; // total number of deaths/births
  Matrix mass;

  particle(){};  
};
