#include "algorithm.hpp"
#include "hamiltonian.hpp"
 

//SEB
void algorithm::Step(particle &X, double dt2){

	Matrix C1;
	C1.set_size(A.row, A.col);
	C1.zeros();

	C1(0, 0) = 1/(1+xi*dt2-A(0, 0)*dt2);
	C1(1, 1) = 1/(1+xi*dt2-A(1, 1)*dt2);
	C1(2, 2) =1/(1+xi*dt2-A(2, 2)*dt2);
	X.p = C1*(X.p + X.f*dt2 + A*X.q*dt2*xi +X.RF3);
	
	X.q = X.q+X.p*dt2; 

	H->compute_force(X);
	H->compute_pressure(X);
	H->add_k_energy(X);

	X.prel = X.p - A*X.q;

}
 