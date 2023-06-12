#ifndef HAMILT_HPP
#define HAMILT_HPP
#include <vector>
#include "matrix.hpp"
#include "particle.hpp"

class hamiltonian
{
public:
	//---------------------------------------------------------------
	//                          Fields
	//---------------------------------------------------------------

	//--- constants for force computations ---
	int np;						// total number of particle = I.Ndim^I.dim
	int dim;					// Number of space dimensions.
	int Ndim_x, Ndim_y, Ndim_z; // particles in each direction.
	int bc;						// boundary condition
	int rep;
	double p1 = 12.;
	double p2 = 6.;		 // Potential of the form r^-p1, r^-p2
	double dcut;	 // cut-off radius for LJ potential
	double sig, eps; // parameters for LJ potential
	double C, D;
	double a;						  // unit cell parameter
	double boxlg_x, boxlg_y, boxlg_z; // periodicity length
	double delx;
	double volume; 
	static constexpr double c = 50.;	 // Strength of the harmonic bonds of the chain.
	static constexpr double kbend = 10.; // Strength of the harmonic bonds of the chain.
										 // Chain potential variables.
	int nChainFreq, nChainLength;

	// simulation box derived parameters. 
	double bottom, top; 
	double mindist, min_box_height;

	//--- Force matrices ---
	Matrix f;
	Matrix f_slow;
	Matrix qtmp; 

	// Used for the KR boundary conditions.
	int reset_count;
	Matrix A, S, Sinv;				  // A (already diagonalized by S)
	Matrix V, Vinv;					  // Eigenvectors for the automorphism group
	Matrix expAt, expmAt;			  // e^At, kept at each step (taking into account BCs)
	Matrix omega, alpha_t, lam, dlam; // We keep track of the exponent in
	double nPer = 0;                      // Remapping number 
	Matrix L, Linv, diag1, offdiag, offdiaginv;

	// r-kr parameters
	double rkr_beta = 1.407715442602953;
	double rkr_mu = 0.281199574322962;
	//--- clist stuff --- 
	int numInteractions;
	int lcmax;
	int nCells;
	double region[3]; // simulation box length in each direction
	int lc[3];		  // number of cells in each direction
	int *head;		  // array representing the set of cells, holding head of linked list
	int *lscl;		  // array, holds id of next particle in same cell at entry

	//---------------------------------------------------------------
	//                        Functions
	//---------------------------------------------------------------

	//--- Kinetic energy ---
	hamiltonian(const input &);
	void add_k_energy(particle &);

	//--- Elementary forces ---
	double LJ(double);
	double dLJ(double);
	double W(double);
	double dW(double);

	//--- global force computations ---
	void compute_force(particle &);
	void compute_force_slow(particle &);
	void compute_pressure(particle &);

	// Functions for the KR boundary conditions.
	void initialize_box(particle &); // Initial positions.
	void Periodic(particle &);
	void deform_box(double);
	double lengthBC(const double, const double, const double,
					double &, double &, double &);

	double CB(double);
	double dCB(double);
	double bend(double);
	double dbend(double);
	double dihedral(double);
	double ddihedral(double);

	// functions for clist
	void initialize_clist();
	void clear_clist();
	void add_particles(particle &);
	void deconstruct_clist();
	void update_force(Matrix &, const double, const double, const double, const double, const int, const int); 

	double chain_length(int, int, particle &);  
	double BCnorm(Matrix); 
	void MinDistance(particle &);
private:
	int index(int, int, int);
	Matrix cross(double, double, double, double, double, double);
	void update_force_r(Matrix &, Matrix, const int);
};

#endif
