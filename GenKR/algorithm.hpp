#ifndef ALGORITH_HPP
#define ALGORITH_HPP

#include "hamiltonian.hpp"
#include "matrix.hpp"
#include <random>

class algorithm
{
public:
  //---------------------------------------------------------------
  //                          Fields
  //---------------------------------------------------------------

  //--- structures ---
  hamiltonian *H;
  particle Xa[5];
  particle Xtmp;

  std::random_device rd;
  int rdd = rd();

  //--- useful constants ---
  double dt, beta, a, xi, sig, dcut;
  int np, Ndim_x, Ndim_y, Ndim_z, dim, rep;
  int n_iter, freq, freq_xmakemol;
  int n_therm_iter;
  double T_therm, T_final;
  Matrix A;     // The strain rate A, which we assume is diagonalized by
                // matrix S.
  Matrix gamma; // Contains twice the symmetrized strain rate tensor
  Matrix qtmp;
  double force_error;
  double numIt;
  double TOL;

  // derived quantities with simulation
  double viscosity, eta_pcf, eta_pef;
  int Nsamp, Nsamp_x, Nsamp_y; // box parameters

  bool genKR;

	algorithm(const input &);
  //---------------------------------------------------------------
  //                        Functions
  //---------------------------------------------------------------

  //---- Initialization ---
  void Initialize();

  //---- Integrators (sampling) ---
  void Step(particle &, double);
  void StepTrue(particle &, double);
  void Therm_Step(particle &, double);

  //--- Main function ---
  void run();

  // Radial Functions
  void RadialObservables(int, ofstream &, Matrix &, int, double);
  void HistogramObservables(ofstream &, Matrix &, int, double, double);

private:
  // Initializations and observables.
  void WritePositions(int, ofstream &);
  void PrepObservablesFile(ofstream &, ofstream &, ofstream &, ofstream &);
  void WriteObservables(double, ofstream &, ofstream &, double, ofstream &, ofstream &, ofstream &);
  int ind(int, int, int = 0);
};
#endif
