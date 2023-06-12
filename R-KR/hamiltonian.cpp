#include "hamiltonian.hpp"
#include "algorithm.hpp"
#include "math.h"

//--- Compute kinetic energy ---
// Add it to the potential energy, which is compted in the force computation..
void hamiltonian::add_k_energy(particle &X)
{
	X.total_momentum_x = 0.;
	X.total_momentum_y = 0.;
	X.total_momentum_z = 0.;
	X.energy = X.potential_energy;
	for (int i = 0; i < np; i++)
	{
		X.total_momentum_x += X.prel(0, i);
		X.total_momentum_y += X.prel(1, i);
		X.total_momentum_z += X.prel(2, i);
		X.energy += 0.5 * X.prel(0, i) * X.prel(0, i) / X.mass(0, i) + 0.5 * X.prel(1, i) * X.prel(1, i) / X.mass(0, i) + 0.5 * X.prel(2, i) * X.prel(2, i) / X.mass(0, i);
	}
}

//------ Regularized Lennard-Jones potential energy --------
double hamiltonian::LJ(double d)
{
	double pot;
	if (d > dcut)
		pot = 0;
	else
		pot = 4. * eps * (pow(sig / d, p1) - pow(sig / d, p2)) + D + C * (d - dcut);
	return pot;
}

//------ Regularized Lennard-Jones force --------
double hamiltonian::dLJ(double d)
{
	double ff;
	if (d > dcut)
		ff = 0;
	else
		ff = -(4. * eps * (p1 * pow(sig / d, p1) - p2 * pow(sig / d, p2)) / d - C);
	return ff;
}

// The chains nearest-neighbor bond potential.
double hamiltonian::CB(double d)
{
	return .5 * c * pow((d - 1), 2);
}

double hamiltonian::dCB(double d)
{
	return c * (d - 1);
}

double hamiltonian::bend(double cos_theta)
{
	double cos_theta0 = -0.36; 
	return .5 * kbend * (cos_theta - cos_theta0) * (cos_theta - cos_theta0);
}
double hamiltonian::dbend(double cos_theta)
{
	double cos_theta0 = -0.36;  
	return kbend * (cos_theta - cos_theta0);
}

double hamiltonian::W(double d)
{
	return pow((d - 1) * (d - 2), 2);
}

double hamiltonian::dW(double d)
{
	return 2 * (d - 1) * (d - 1) * (d - 2) + 2 * (d - 1) * (d - 2) * (d - 2);
}

void hamiltonian::compute_pressure(particle &X)
{
	X.pressure.zeros();
	for (int i = 0; i < np; i++)
	{
		// Always use 3 dimensions, even if z is always zero.
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				X.pressure(j, k) += X.prel(j, i) * X.prel(k, i) / X.mass(0, i);
			}
		}
	}
	X.pressure = X.pressure * (1. / volume);
}

void hamiltonian::update_force_r(Matrix &f, Matrix g, const int i)
{
	f(0, i) += g(0, 0);
	f(1, i) += g(1, 0);
	if (Ndim_z > 1)
	{
		f(2, i) += g(2, 0);
	}
}

void hamiltonian::update_force(Matrix &f, const double act, const double dx, const double dy, const double dz, const int i, const int j)
{
	f(0, i) += act * dx;
	f(1, i) += act * dy;
	f(0, j) -= act * dx;
	f(1, j) -= act * dy;
	if (Ndim_z > 1)
	{
		f(2, i) += act * dz;
		f(2, j) -= act * dz;
	}
}

// force on all the particles, plus the potential energy of the system
void hamiltonian::compute_force(particle &X)
{ 
	double act;
	double dist;
	double dx, dy, dz;
	numInteractions = 0;

	int observHold = 0;

	//--- reset forces to 0 ---
	X.f.zeros();
	X.sigma.zeros();
	double enpot = 0.;

	clear_clist();
	add_particles(X);

	double dist2, dist3;
	double dx2, dy2, dz2, dx3, dy3, dz3;
	// Chain interactions.
	for (int i = 0; i < np; i++)
	{
		if ((i % nChainFreq) < nChainLength - 1 && i + 1 < np)
		{
			dist = lengthBC(X.q(0, i) - X.q(0, i + 1), X.q(1, i) - X.q(1, i + 1),
							X.q(2, i) - X.q(2, i + 1), dx, dy, dz);
			if (dist > min_box_height / 2)
			{
				cout << "ERROR: chain bond length longer than box height/2." << endl
					 << "box height " << min_box_height << " dist " << dist
					 << " particles " << i << " and " << i + 1 << endl;
				exit(-1);
			}
			enpot += CB(dist);
			act = -dCB(dist);
			update_force(X.f, act, dx, dy, dz, i, i + 1);
		}
		if ((i % nChainFreq) < nChainLength - 2 && i + 2 < np)
		{
			// Note, we already have dist from previous if statement.
			dist2 = lengthBC(X.q(0, i + 1) - X.q(0, i + 2), X.q(1, i + 1) - X.q(1, i + 2),
							 X.q(2, i + 1) - X.q(2, i + 2), dx2, dy2, dz2);
			// Double bond energy
			double cos_theta = -(dx * dx2 + dy * dy2 + dz * dz2);
			enpot += bend(cos_theta);
			act = -dbend(cos_theta);
			update_force(X.f, -act / dist, dx2, dy2, dz2, i, i + 1);
			update_force(X.f, -act * cos_theta / dist, dx, dy, dz, i, i + 1);
			update_force(X.f, -act / dist2, dx, dy, dz, i + 1, i + 2);
			update_force(X.f, -act * cos_theta / dist2, dx2, dy2, dz2, i + 1, i + 2);
		}
	}
	//  Build cell list

	int c;		// for storing cell index where particle resides
	int adjC;	// for storing cell index of adjacent cell
	int mc[3];	// holds x,y,z index of a cell
	int mcl[3]; // a handy mc[a] for adjacent cells
	double da[3];  
	int j; // for going through linked lists
	int testHold = 0;

	// Solvent-solvent interactions.
	for (int i = 0; i < np; i++)
	{
		da[0] = Linv(0, 0) * X.q(0, i) + Linv(0, 1) * X.q(1, i) + Linv(0, 2) * X.q(2, i);
		da[1] = Linv(1, 0) * X.q(0, i) + Linv(1, 1) * X.q(1, i) + Linv(1, 2) * X.q(2, i);
		da[2] = Linv(2, 0) * X.q(0, i) + Linv(2, 1) * X.q(1, i) + Linv(2, 2) * X.q(2, i);

		for (int a = 0; a < 3; a++)
		{
			da[a] += .5;
		} 
		// getting vector cell indexes
		for (int a = 0; a < 3; a++)
		{
			mc[a] = (int)floor(abs(da[a] * lc[a]));
			mc[a] = (mc[a] + lc[a]) % lc[a];
		}
		c = mc[0] * lc[1] * lc[2] + mc[1] * lc[2] + mc[2];

		// go through adjacent cells
		for (mcl[0] = mc[0] - 1; mcl[0] <= mc[0] + 1; mcl[0]++)
		{
			for (mcl[1] = mc[1] - 1; mcl[1] <= mc[1] + 1; mcl[1]++)
			{
				// Special handling in case lc[2] = 1, 2 (to handle planar flow).
				int zcm = (lc[2] > 2 ? 1 : 0);
				int zcp = (lc[2] > 1 ? 1 : 0);
				for (mcl[2] = mc[2] - zcm; mcl[2] <= mc[2] + zcp; mcl[2]++)
				{
					adjC = ((mcl[0] + lc[0]) % lc[0]) * lc[1] * lc[2] + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[2] + lc[2]) % lc[2]);

					if (adjC < 0 || adjC >= lc[0] * lc[1] * lc[2])
					{
						cout << "ERROR: adjC out of bounds " << adjC << endl;
						exit(-1);
					}

					j = head[adjC];
					while (j != -1)
					{
						// avoid double counting pairs
						if (i < j)
						{
							// compute forces
							numInteractions++;
							dist = lengthBC(X.q(0, i) - X.q(0, j), X.q(1, i) - X.q(1, j), X.q(2, i) - X.q(2, j), dx, dy, dz);
							enpot += LJ(dist);
							act = -dLJ(dist);
							update_force(X.f, act, dx, dy, dz, i, j);

							mindist = min(dist, mindist);
							X.sigma(0, 0) += (dist * dx) * (act * dx);
							X.sigma(1, 0) += (dist * dy) * (act * dx);
							X.sigma(2, 0) += (dist * dz) * (act * dx);
							X.sigma(1, 1) += (dist * dy) * (act * dy);
							X.sigma(2, 1) += (dist * dz) * (act * dy);
							X.sigma(2, 2) += (dist * dz) * (act * dz);
						}
						j = lscl[j];
					}
				}
			}
		}
	}
	X.sigma(0, 1) = X.sigma(1, 0);
	X.sigma(0, 2) = X.sigma(2, 0);
	X.sigma(1, 2) = X.sigma(2, 1);

	// potential energy of replica 'k'
	X.energy = enpot;
	X.potential_energy = enpot;

	X.sigma = X.sigma * (1. / volume);
}
//--------- end of force computations ----









//-------- Apply Generalized KR boundary conditions ---------
void hamiltonian::Periodic(particle &X)
{
	Matrix qtmp = Linv * X.q;
	double ipart;

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < np; j++)
		{
			qtmp(i, j) = qtmp(i, j) - round(qtmp(i, j));
		}
	}
	X.q = L * qtmp;
	X.p = X.prel + A * X.q; // Recover absolute velocity
} 


//-------- Minimum Distance ---------
void hamiltonian::MinDistance(particle &X)
{
	double dmin = 1000.0;
	int minL = -4;
	int maxL = 4;
	double dx, dy, dz, dist;
	for (int i = minL; i < maxL; i++)
	{
		for (int j = minL; j < maxL; j++)
		{
			for (int k = minL; k < maxL; k++)
			{ 
				if (i != 0 || j != 0 || k != 0)
				{
					dx = L(0, 0) * i + L(0, 1) * j + L(0, 2) * k;
					dy = L(1, 0) * i + L(1, 1) * j + L(1, 2) * k;
					dz = L(2, 0) * i + L(2, 1) * j + L(2, 2) * k;
					dist = sqrt(dx * dx + dy * dy + dz * dz + 1e-32); // distance up to BC
					if (dist < dmin)
					{
						dmin = dist;
					}
				}
			}
		}
	}
	X.mindist = dmin ;
}


///-----------------------------------------------
//          Elementary interactions
//------------------------------------------------
// This version is the unrolled loop version
double hamiltonian::lengthBC(
	const double dxin, const double dyin, const double dzin,
	double &dx, double &dy, double &dz)
{
	double dist;
	// We transfer into lattice basis coordinates.  Note that
	// this computes the minimum periodic distance only if
	// that distance is less than half the smallest height of the
	// lattice cell.
	//
	// If rcut is less than 1/2 *min_i h_i for all time (reasonable
	// due to boundedness of the deformation), then the error we
	// are making is not a problem!
	dx = Linv(0, 0) * dxin + Linv(0, 1) * dyin + Linv(0, 2) * dzin;
	dy = Linv(1, 0) * dxin + Linv(1, 1) * dyin + Linv(1, 2) * dzin;
	dz = Linv(2, 0) * dxin + Linv(2, 1) * dyin + Linv(2, 2) * dzin;
	// Reference coordinates
	double tdx = dx - round(dx);
	double tdy = dy - round(dy);
	double tdz = dz - round(dz);

	dx = L(0, 0) * tdx + L(0, 1) * tdy + L(0, 2) * tdz;
	dy = L(1, 0) * tdx + L(1, 1) * tdy + L(1, 2) * tdz;
	dz = L(2, 0) * tdx + L(2, 1) * tdy + L(2, 2) * tdz;
	dist = sqrt(dx * dx + dy * dy + dz * dz + 1e-32); // distance up to BC
	dx = dx / dist;
	dy = dy / dist;
	dz = dz / dist;
	return dist;
}

double hamiltonian::BCnorm(Matrix m)
{ // only for positions matrix
	Matrix dist;
	dist.set_size(1, m.col);
	dist.zeros();
	double dx, dy, dz;
	double dxin, dyin, dzin;
	for (int i = 0; i < m.col; i++)
	{
		dx = m(0, i);
		dy = m(1, i);
		dz = m(2, i);
		dist(0, i) = lengthBC(dx, dy, dz, dxin, dyin, dzin);
	}

	return dist.norm();
}


double hamiltonian::chain_length(int start, int end, particle &X)
{
	double l_temp, qx, qy, qz;
	double dx, dy, dz;
	qx = 0;
	qy = 0;
	qz = 0;
	for (int i = start; i < end; i++)
	{
		l_temp = lengthBC(X.q(0, i) - X.q(0, i + 1), X.q(1, i) - X.q(1, i + 1), X.q(2, i) - X.q(2, i + 1), dx, dy, dz);
		qx += l_temp * dx;
		qy += l_temp * dy;
		qz += l_temp * dz;
	}
	return sqrt(qx * qx + qy * qy + qz * qz);
}



void hamiltonian::initialize_box(particle &X)
{
	int Ndim_x_new = Ndim_x;
	int Ndim_y_new = Ndim_y;
	Matrix qlam(3, 1); // New particle position in lattice coordinates.

	for (int l = 0; l < Ndim_z; ++l)
	{
		for (int i = 0; i < Ndim_x_new; ++i)
		{
			for (int j = 0; j < Ndim_y_new; ++j)
			{
				if (i * Ndim_y_new + j < Ndim_x * Ndim_y)
				{
					qlam(0, 0) = (0.5 + i - 0.5 * Ndim_x_new) / Ndim_x_new;
					qlam(1, 0) = (0.5 + j - 0.5 * Ndim_y_new) / Ndim_y_new;
					qlam(2, 0) = (0.5 + l - 0.5 * Ndim_z) / Ndim_z;
					int ind = (i * Ndim_y_new + j) * Ndim_z + l;
					qlam = L * qlam;
					X.q(0, ind) = qlam(0, 0);
					X.q(1, ind) = qlam(1, 0);
					X.q(2, ind) = qlam(2, 0);
				}
			}
		}
	}
}






