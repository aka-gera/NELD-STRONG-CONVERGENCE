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
 