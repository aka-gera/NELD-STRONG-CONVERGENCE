#include "algorithm.hpp"
#include "hamiltonian.hpp" 

 
  //--- Creation: initialize some fields ---
  algorithm::algorithm(const input &I)
  {
    rep = 5;
    TOL = 1e-9;
    np = I.Ndim_x * I.Ndim_y * I.Ndim_z;
    beta = I.beta;
    a = I.a;
    Ndim_x = I.Ndim_x;
    Ndim_y = I.Ndim_y;
    Ndim_z = I.Ndim_z;
    dim = I.dim;
    dt = I.t_step;
    T_therm = I.T_therm;
    T_final = I.T_final;
    xi = I.xi;
    n_iter = static_cast<int>(std::ceil(I.T_final / dt));                           
    n_therm_iter = static_cast<int>(std::ceil(I.T_therm / dt));                     
    freq = static_cast<int>(std::max(static_cast<double>(n_iter) / I.n_obs, 1.0));  
    if (I.n_xmakemol > 0)
    {
      freq_xmakemol = static_cast<int>(std::max(static_cast<double>(n_iter) / I.n_xmakemol, 1.0));  
    }
    else
    {
      freq_xmakemol = -1;
    }
    sig = I.sig;
    dcut = I.dcut;
    A = I.A;

    // Following Baranyai and Cummings, I do not divide by 2.
    gamma = A + A.transpose();
    genKR = I.genKR;
    Nsamp_x = 20;
    Nsamp_y = 20;
    Nsamp = Nsamp_x * Nsamp_y; // Just record the quadrants.
  }


void algorithm::Initialize()
{
	std::mt19937 gen(rdd);
	std::uniform_real_distribution<> rand_u(0., 1.);
	std::normal_distribution<> rand_n(0., 1.);

	int seed = rdd;
	qtmp.set_size(dim, np);
	qtmp.zeros();

	//--- Initialize the fields ---
	Xtmp.q.set_size(dim, np);
	Xtmp.q.zeros();
	Xtmp.mass.set_size(1, np);
	Xtmp.mass.zeros();
	Xtmp.p.set_size(dim, np);
	Xtmp.p.zeros();
	Xtmp.RF.set_size(dim, np);
	Xtmp.RF.zeros();
	Xtmp.f.set_size(dim, np);
	Xtmp.f.zeros();
	Xtmp.pressure.set_size(dim, dim);
	Xtmp.pressure.zeros();
	Xtmp.sigma.set_size(dim, dim);
	Xtmp.sigma.zeros();
	Xtmp.energy = 0;
	Xtmp.potential_energy = 0;

	cout << "Starting with seed " << seed << endl;

	for (int r = 0; r < rep; ++r)
	{
		//--- Initialize the fields ---
		Xa[r].q.set_size(dim, np);
		Xa[r].q.zeros();
		Xa[r].mass.set_size(1, np);
		Xa[r].mass.zeros();
		Xa[r].p.set_size(dim, np);
		Xa[r].p.zeros();
		Xa[r].RF.set_size(dim, np);
		Xa[r].RF.zeros();
		Xa[r].RF2.set_size(dim, np);
		Xa[r].RF2.zeros();
		Xa[r].RF3.set_size(dim, np);
		Xa[r].RF3.zeros();
		Xa[r].RF4.set_size(dim, np);
		Xa[r].RF4.zeros();
		Xa[r].f.set_size(dim, np);
		Xa[r].f.zeros();
		Xa[r].pressure.set_size(dim, dim);
		Xa[r].pressure.zeros();
		Xa[r].sigma.set_size(dim, dim);
		Xa[r].sigma.zeros();
		Xa[r].energy = 0;
		Xa[r].potential_energy = 0;

		H->initialize_box(Xa[r]);

		Xa[r].p = A * Xa[r].q;
	}

	int effdim = (Ndim_z > 1 ? 3 : 2);
	for (int i = 0; i < effdim; ++i)
	{
		for (int j = 0; j < np; ++j)
		{
			Xa[0].mass(0, j) = 1.0;
			Xa[0].q(i, j) += 0.05 * rand_n(gen);
			Xa[0].p(i, j) += sqrt(1. / beta) * rand_n(gen);
			for (int r = 1; r < rep; ++r)
			{
				Xa[r].mass(0, j) = 1.0;
				Xa[r].q(i, j) = Xa[0].q(i, j);
				Xa[r].p(i, j) = Xa[0].p(i, j);
			}
		}
	}

	for (int r = 0; r < rep; ++r)
	{
		Xa[r].prel = Xa[r].p - A * Xa[r].q;
		// Normalize to the correct initial temperature.
		double alpha = sqrt(1. / beta * dim * np / Xa[r].prel.contraction_product(Xa[r].prel));
		Xa[r].p -= Xa[r].prel * (1. - alpha);
		Xa[r].prel = Xa[r].p - A * Xa[r].q;
		H->Periodic(Xa[r]);
	}

	ofstream POUT("out_param.dat");
	//------------------- Screen output --------------------------
	POUT << endl;
	POUT << "------------------------------------------" << endl;
	POUT << "     CONSTITUTIVE LAW COMPUTATION         " << endl;
	POUT << "------------------------------------------" << endl;
	POUT << endl;
	POUT << "        Initialization performed          " << endl;
	POUT << "                                          " << endl;
	POUT << "-----  Parameters of the computation -----" << endl;
	POUT << endl;
	POUT << " Time step                 (dt) : " << dt << endl;
	POUT << " Number of iterations           : " << n_iter << endl;
	POUT << " Xmakemol output every...steps  : " << freq_xmakemol << endl;
	POUT << " Other outputs   every...steps  : " << freq << endl;
	POUT << " Friction coefficient      (xi) : " << xi << endl;
	POUT << " Inverse temperature     (beta) : " << beta << endl;
	POUT << " Temperature           (1/beta) : " << 1 / beta << endl;
	POUT << " Spatial dimension        (dim) : " << dim << endl;
	POUT << " Total number of particles (np) : " << np << endl;
	POUT << " Elementary cell size       (a) : " << a << endl;
	POUT << " Particle density               : " << 1 / (a * a * a) << endl;
	POUT << " WCA equilibrium distance (sig) : " << sig << endl;
	POUT << " WCA energy epsilon       (eps) : " << H->eps << endl;
	POUT << " Cut off radius   (sigma units) : " << dcut << endl;
	POUT << " Background flow            (A) : " << A << endl;
	POUT << " Random number seed      (seed) : " << seed << endl;

	POUT << endl;
	POUT << "----------------------------------------- " << endl;

	POUT << endl;
}



// SOILE-B 

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

	Gamma(0, 0) = xi - A(0, 0);
	Gamma(1, 1) = xi - A(1, 1);
	Gamma(2, 2) = xi - A(2, 2);

	X.p = X.p + (X.f + A*X.q*xi - Gamma*X.p)*dt*0.5 - Gamma*((X.f + A*X.q*xi - Gamma*X.p)*dt*dt*0.5 + (X.RF2 - X.RF*dt*0.25) * 2)*0.25 + X.RF*0.5;

	X.q = X.q + X.p*dt + X.RF2 - X.RF*dt*0.5;

	X.prel = X.p - A*X.q;
	Xtmp = X;

	H->Periodic(Xtmp);
	H->compute_force(Xtmp);
	X.f = Xtmp.f;

	H->compute_pressure(X);
	H->add_k_energy(X);

	X.p = X.p + (X.f + A*X.q*xi - Gamma*X.p)*dt*0.5 - Gamma*((X.f + A*X.q*xi - Gamma*X.p)*dt*dt*0.5 + (X.RF2 - X.RF*dt*0.25) * 2)*0.25 + X.RF*0.5;

	X.prel = X.p - A*X.q;
} 




void algorithm::Therm_Step(particle &X, double dt)
{
	X.p = (X.p + X.f * dt) * (1 / (1 + xi * dt));
	X.prel = X.p - A * X.q;
	X.q += X.prel * dt; // TODO: only works for mass==1 ! 
	X.p += X.RF + X.RF2 + X.RF3;  
	X.prel = X.p - A * X.q;

	H->Periodic(X);
	H->compute_force(X); 
	H->compute_pressure(X);
	H->add_k_energy(X);
}
 


void algorithm::WritePositions(int Niter, ofstream &XMAKEMOL)
{
	int count = max(rep * H->nChainLength, rep);
	char types[] = "HCNOPS";
	XMAKEMOL << count << endl;
	XMAKEMOL << "Time = " << Niter * dt << endl;

	for (int r = 0; r < rep; r++)
	{
		XMAKEMOL << types[r % 6] << " ";
		double qx = Xa[0].q(0, 0);
		double qy = Xa[0].q(1, 0);
		double qz = Xa[0].q(2, 0);
		double dx, dy, dz, dist;
		if (r > 0)
		{
			dist = H->lengthBC(Xa[r].q(0, 0) - Xa[0].q(0, 0), Xa[r].q(1, 0) - Xa[0].q(1, 0), Xa[r].q(2, 0) - Xa[0].q(2, 0), dx, dy, dz);
			qx += dist * dx;
			qy += dist * dy;
			qz += dist * dz;
		}
		XMAKEMOL << qx << "  " << qy << "  " << qz << "  ";
		for (int i = 1; i < H->nChainLength; ++i)
		{
			dist = H->lengthBC(Xa[r].q(0, i) - Xa[r].q(0, i - 1), Xa[r].q(1, i) - Xa[r].q(1, i - 1), Xa[r].q(2, i) - Xa[r].q(2, i - 1), dx, dy, dz);
			XMAKEMOL << " " << dx * dist << " " << dy * dist << " " << dz * dist << endl;

			XMAKEMOL << types[r % 6] << " ";
			qx += dist * dx;
			qy += dist * dy;
			qz += dist * dz;
			XMAKEMOL << qx << "  " << qy << "  " << qz << "  ";
		}
		XMAKEMOL << endl;
	}

} // WritePositions()

void algorithm::PrepObservablesFile(ofstream &ENERGY, ofstream &ENERGYDIFF, ofstream &QNORMS, ofstream &PNORMS)
{
	ENERGY << "#" << setw(11) << "time"
		   << " ";
	for (int r = 0; r < rep; r++)
	{
		ENERGY << setw(12) << "Potential"
			   << " "
			   << setw(12) << "Kinetic"
			   << " ";
	}
	ENERGY << endl;
	ENERGYDIFF << "#" << setw(11) << "time"
			   << " ";
	for (int r = 0; r < rep; r++)
	{
		ENERGYDIFF << setw(12) << "Potential"
				   << " "
				   << setw(12) << "Kinetic"
				   << " ";
	}
	ENERGYDIFF << endl;
	QNORMS << "#" << setw(11) << "time"
		   << " ";
	QNORMS << setw(20) << "Xa[0].q - Xa[0].q"
		   << " "; 
	for (int r = 1; r < rep; r++)
	{
		QNORMS << setw(6) << "Xa[" << r << "].q - Xa[" << 0 << "].q"
			   << " ";
	}
	QNORMS << endl;
	PNORMS << "#" << setw(11) << "time"
		   << " ";
	PNORMS << setw(20) << "Xa[0].p - Xa[0].p"
		   << " "; 
	for (int r = 1; r < rep; r++)
	{
		PNORMS << setw(6) << "Xa[" << r << "].p - Xa[" << 0 << "].p"
			   << " ";
	}
	PNORMS << endl;
}

void algorithm::WriteObservables(double time, ofstream &ENERGY, ofstream &ENERGYDIFF, double mean_kin, ofstream &QNORMS, ofstream &PNORMS, ofstream &MINDIST)
{
	ENERGY << setw(12) << time << " ";
	for (int r = 0; r < rep; r++)
	{
		ENERGY << setw(12) << (Xa[r].potential_energy) << " "
			   << setw(12) << (Xa[r].energy - Xa[r].potential_energy) << " ";
	}
	ENERGY << endl;

	MINDIST << setw(12) << time << " ";
	for (int r = 0; r < rep; r++)
	{
		MINDIST << setw(12) << Xa[r].mindist << " ";
	}
	MINDIST << endl;

	ENERGYDIFF << setw(12) << time << " ";

	for (int r = 0; r < rep; r++)
	{
		ENERGYDIFF << setw(12) << (Xa[0].potential_energy - Xa[r].potential_energy) << " "
				   << setw(12) << ((Xa[0].energy - Xa[0].potential_energy) - (Xa[r].energy - Xa[r].potential_energy)) << " ";
	}

	ENERGYDIFF << endl;

	QNORMS << setw(12) << time << " " << setw(20) << (Xa[0].q - Xa[0].q).norm() << " ";

	for (int r = 1; r < rep; r++)
	{
		double tdx, tdy, tdz;
		int ret = 0;

		QNORMS << setw(25) << H->BCnorm(Xa[r].q - Xa[0].q) << " ";
	}
	QNORMS << endl;

	PNORMS << setw(12) << time << " " << setw(20) << (Xa[0].p - Xa[0].p).norm() << " ";

	for (int r = 1; r < rep; r++)
	{
		PNORMS << setw(25) << (Xa[r].prel - Xa[0].prel).norm() << " ";
	}
	PNORMS << endl;
}
 