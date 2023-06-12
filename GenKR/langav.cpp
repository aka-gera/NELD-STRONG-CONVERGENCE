#include "algorithm.hpp"
#include "hamiltonian.hpp"

// EE
void algorithm::Step(particle &X, double dt2)
{
	qtmp = X.q + X.p * dt2;
	H->compute_force(X);
	H->compute_pressure(X);
	H->add_k_energy(X);
	X.p += X.f * dt2 - X.p * (dt2 * xi) + A * X.q * dt2 * xi + A * X.p * dt2 + X.RF3;
	X.q = qtmp;
	X.prel = X.p - A * X.q;
}

// SOILE-B
void algorithm::StepTrue(particle &X, double dt)
{

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

	C2 = (X.f + A * X.q * xi - Gamma * X.p) * dt * dt * 0.5 + X.RF2; 
	Ftmp = (X.f + A * X.q * xi); 
	X.q = X.q + X.p * dt + C2;

	X.prel = X.p - A * X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);
	X.p = X.p + (Ftmp + X.f + A * X.q * xi) * dt * 0.5 - Gamma * (X.p * dt + C2) + X.RF;

	X.prel = X.p - A * X.q;
}
