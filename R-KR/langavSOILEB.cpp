#include "algorithm.hpp"
#include "hamiltonian.hpp"




//SOILE-B
void algorithm::Step(particle &X, double dt2){

	Matrix Gamma2;
	Gamma2.set_size(A.row, A.col);
	Gamma2.zeros();
 

	Matrix Ftmp;
	Ftmp.set_size(A.row, 1);
	Ftmp.zeros();

	Gamma2(0,0) = xi - A(0,0);
	Gamma2(1,1) = xi - A(1,1);
	Gamma2(2,2) = xi - A(2,2);

	X.p = X.p + (X.f + A*X.q*xi - Gamma2*X.p)*dt2*0.5 - Gamma2*((X.f + A*X.q*xi - Gamma2*X.p)*dt2*dt2*0.5 + (X.RF4 - X.RF3*dt2*0.25)*2)*0.25 + X.RF3*0.5;

	X.q = X.q + X.p*dt2 +X.RF4-X.RF3*dt2*0.5;

	X.prel = X.p - A*X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);

	X.p = X.p + (X.f + A*X.q*xi - Gamma2*X.p)*dt2*0.5 - Gamma2*((X.f + A*X.q*xi - Gamma2*X.p)*dt2*dt2*0.5  + (X.RF4 - X.RF3*dt2*0.25)*2)*0.25 + X.RF3*0.5;

	X.prel = X.p - A*X.q;

}


