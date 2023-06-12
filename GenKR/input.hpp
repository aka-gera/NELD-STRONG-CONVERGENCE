//-----------------------------------------------------
//
//     Reading parameters from "input_file_..."
//
//-----------------------------------------------------


#ifndef INPUT_HPP
#define INPUT_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

#include "matrix.hpp"
//--------------- definition of the structure --------
class input
{
	public: 
		//------------------------ FIELDS --------------------------

		//--- Simulation parameters --- 
		double t_step;         // time step for the time integrator algorithm
		double T_therm;
		double T_final; 
		int    n_xmakemol;  // XMakeMol output frequency
		int    n_obs;           // outputs written every 'freq' steps (energy, reaction coordinate, ...)
		double xi;             // friction coefficient for Langevin dynamics
		double beta;        // inverse temperature 
		double  st_rt_a,st_rt_b,st_rt_c;

		//--- Thermodynamic conditions ---
		int dim;  //Spatial particles.
		int space_dim;
		int    Ndim_x,Ndim_y,Ndim_z;        // Number of particles in each direction
		double a;           // domain [0,Ndim*a]^dim, volume = (Ndim*a)^dim 
		bool genKR;  
		Matrix A, S, Sinv;
		//--- Potential parameters ---
		double sig, eps,dcut;   
		//---------------------- READING FUNCTIONS --------------------------

		template <class T>
		//--- elementary reading function -------
		void read_item(istream& ,const char *, T *);
		//--- copying datas from input_file to the parameters class ---
		void load(void);  
};

#endif


