#include "algorithm.hpp"
#include "hamiltonian.hpp"
 

//SEAC
void algorithm::Step(particle &X, double dt2){

	X.q = X.q + X.p*dt2;
	X.prel = X.p - A*X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);

	X.p = X.p + X.f*dt2 - X.p*(dt2*xi) + A*X.q*dt2*xi + A*X.p*dt2 + X.RF3;

	X.prel = X.p - A*X.q; 

}



