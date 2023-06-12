#include "algorithm.hpp"
#include "hamiltonian.hpp"







//ABAO
void algorithm::Step(particle &X, double dt2){

	Matrix C1;
	C1.set_size(A.row, A.col);
	C1.zeros();

	C1(0, 0) = exp(A(0, 0)*dt2);
	C1(1, 1) = exp(A(1, 1)*dt2);
	C1(2, 2) = exp(A(2, 2)*dt2);


	X.p += X.f*0.5*dt2;
	X.q += X.p*dt2; //TODO: only works for mass==1 !
	X.prel = X.p - A*X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);

	double decay = exp(-xi*dt2);
	double diffuse = sqrt((1 - decay*decay) / beta);

	X.p += X.f*0.5*dt2;
	//	X.p = (A*dt).exp_diag()*X.p;
	X.p = C1*X.p;

	X.p = X.p*decay + A*X.q*(1 - decay) + X.RF3;
	X.prel = X.p - A*X.q;

	H->compute_pressure(X);
	//--- Add kinetic energy to the total energy ---
	H->add_k_energy(X);

}





//SOILE-B
void algorithm::StepTrue(particle &X, double dt){

	Matrix Gamma;
	Gamma.set_size(A.row, A.col);
	Gamma.zeros();

	Matrix C2;
	C2.set_size(A.row, 1);
	C2.zeros();

	Matrix Ftmp;
	Ftmp.set_size(A.row, 1);
	Ftmp.zeros();

	Gamma(0, 0) = (xi - A(0, 0));
	Gamma(1, 1) = (xi - A(1, 1));
	Gamma(2, 2) = (xi - A(2, 2));


	C2 = (X.f + A*X.q*xi - Gamma*X.p)*dt*dt*0.5 + X.RF2;



	Ftmp = (X.f + A*X.q*xi);

	X.q = X.q + X.p*dt + C2;

	X.prel = X.p - A*X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);
	X.p = X.p + (Ftmp + X.f + A*X.q*xi)*dt*0.5 - Gamma*(X.p*dt + C2) + X.RF;


	X.prel = X.p - A*X.q;
}


