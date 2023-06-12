#include "algorithm.hpp"
#include "hamiltonian.hpp"
#include <time.h>
#include <iostream>

//----------------------------------------------------------
//     Sampling procedure = compute the trajectory
//----------------------------------------------------------


void algorithm::run()
{
	std::mt19937 gen(rdd);
	std::uniform_real_distribution<> rand_u(0., 1.);
	std::normal_distribution<> rand_n(0., 1.);

	//--------------- Useful fields ---------------------
	ofstream ENERGY("energy.dat");			// observables file
	ofstream XMAKEMOL("xmakemol.xyz");		// XMakeMol output
	ofstream QNORMS("ql2norms.dat");		// l2 norms output of X.q
	ofstream PNORMS("pl2norms.dat");		// l2 norms output of X.p
	ofstream ENERGYDIFF("energy_diff.dat"); // energy differences from finest stepsize
	ofstream MINDIST("mindist.dat");

	//------------- Initialization ----------------------
	for (int r = 0; r < rep; ++r)
	{
		H->compute_force(Xa[r]);
		H->add_k_energy(Xa[r]);
	}

	int acquisition = 0;
	int acquisition_xmakemol = 0;
	int smfreq = max(1, freq / 10);

	cout << "------- Computation started -------" << endl;
	cout << endl;
	if (freq_xmakemol > 0)
	{
		WritePositions(0, XMAKEMOL);
	}
	PrepObservablesFile(ENERGY, ENERGYDIFF, QNORMS, PNORMS);

	double kinetic;
	double mean_kinetic = 0;
	double fluct = sqrt(2 * dt * xi / beta);
	double fluct2 = sqrt(2 * dt * dt * dt * xi / beta);
	int effdim = (Ndim_z > 1 ? 3 : 2);

	double tdt = 0.001;
	double decay = exp(-xi * dt);
	double diffuse = sqrt((1 - decay * decay) / beta);
	// Thermalization to equilibrium.  Here we take rather large steps
	n_therm_iter = 1000;

	for (int Niter = -n_therm_iter; Niter < 0; Niter++)
	{
		for (int i = 0; i < np; ++i)
		{
			for (int j = 0; j < effdim; ++j)
			{
				Xa[0].RF(j, i) = rand_n(gen) * fluct;
			}
		}
		Therm_Step(Xa[0], tdt);

		Xa[0].RF.zeros();
		Xa[0].RF2.zeros();
		Xa[0].RF3.zeros();
		Xa[0].RF4.zeros();
		if (Niter % ((n_therm_iter) / 25) == 0)
		{
			WriteObservables(Niter * tdt, ENERGY, ENERGYDIFF, mean_kinetic, QNORMS, PNORMS, MINDIST);
			cout << 100 + Niter * 100 / n_therm_iter << "% done.  "
				 << setw(12) << Xa[0].potential_energy << " "
				 << setw(12) << Xa[0].energy - Xa[0].potential_energy << endl;
		}
	}
	Xa[0].p = Xa[0].prel + A * Xa[0].q;
	for (int r = 0; r < rep; ++r)
	{
		Xa[r].q = Xa[0].q;
		Xa[r].p = Xa[0].p;
		Xa[r].prel = Xa[0].prel;
		Xa[r].mindist = Xa[0].mindist;
		H->compute_force(Xa[r]);
		H->compute_pressure(Xa[r]);
		H->add_k_energy(Xa[r]);
	}

	WriteObservables(0., ENERGY, ENERGYDIFF, mean_kinetic, QNORMS, PNORMS, MINDIST);
	//-------------- Loop over time ---------------------
	cout << "Main Simulation " << endl;
	for (int Niter = 1; Niter < n_iter + 1; Niter++)
	{
		for (int i = 0; i < np; ++i)
		{
			for (int j = 0; j < effdim; ++j)
			{
				double xx = rand_n(gen);
				double yy = rand_n(gen);
				Xa[0].RF(j, i) = xx * fluct;
				Xa[0].RF2(j, i) = (xx + yy * (1 / sqrt(3))) * fluct2 * 0.5;
				Xa[0].RF3(j, i) = xx * fluct;
				Xa[0].RF4(j, i) = (xx + yy * (1 / sqrt(3))) * fluct2 * 0.5;
				for (int k = 1; k < rep; ++k)
				{
					Xa[k].RF4(j, i) = Xa[0].RF4(j, i) + Xa[k].RF4(j, i) + Xa[k].RF3(j, i) * dt;
					Xa[k].RF3(j, i) = Xa[0].RF3(j, i) + Xa[k].RF3(j, i);
				}
			}
		}

		StepTrue(Xa[0], dt);
		H->Periodic(Xa[0]);
		H->MinDistance(Xa[0]);
		Xa[0].RF.zeros();
		Xa[0].RF2.zeros();

		for (int r = 1; r < rep; ++r)
		{
			if (Niter % ((int)pow(2, r)) == 0)
			{ 
				Step(Xa[r], pow(2, r) * dt);
				H->Periodic(Xa[r]);
				Xa[r].RF3.zeros();
				Xa[r].RF4.zeros();
			}
		}

		H->deform_box(dt);
		acquisition_xmakemol += 1;
		if (Niter % ((n_iter + 1) / 10) == 0)
		{
			cout << Niter * 10 / n_iter << "0% done.  "
				 << setw(12) << Xa[0].potential_energy << " "
				 << setw(12) << Xa[0].energy - Xa[0].potential_energy << endl;
		}

		if (Niter % freq == 0)
		{
			WriteObservables(Niter * dt, ENERGY, ENERGYDIFF, mean_kinetic, QNORMS, PNORMS, MINDIST);
			for (int r = 0; r < rep; ++r)
			{
				Xa[r].mindist = Xa[0].mindist;
			}
		}

		//---- Positions (read with XMakeMol software) ------
		if (acquisition_xmakemol == freq_xmakemol)
		{
			acquisition_xmakemol = 0;
			WritePositions(Niter, XMAKEMOL);
		}
	}
	cout << endl;
	for (int r = 0; r < rep; r++)
	{
		cout << "(Xa[0].q - Xa[" << r << "].q).norm() = " << (Xa[0].q - Xa[r].q).norm() << endl;
		cout << "(Xa[0].p - Xa[" << r << "].p).norm() = " << (Xa[0].p - Xa[r].p).norm() << endl;
		cout << endl;
	}

	//------------- End of computation ---------------
	cout << endl;
	cout << "------- End of computation ------" << endl;
	cout << endl;
	cout << endl;
	cout << "Output files:" << endl;
	cout << "  energy.dat:  time step, total energy, potential energy, "
		 << "avg kinetic" << endl;
	cout << "  ql2norms.dat and pl2norms.dat gives norm outputs of X.q and X.p" << endl;
	cout << "  radial.dat " << endl;
	cout << "  torsion.dat " << endl;
	cout << "  xmakemol.xyz:  positions of full system" << endl;
	cout << "  length.dat gives a histogram of the distance between the first and last particle in the chain" << endl;
	cout << "  dimer.dat: outputs each step of the bond length" << endl;
	cout << " BG Flow " << A << endl;
	cout << "Total reset count for the simulation " << H->reset_count << endl;
	cout << endl;

	QNORMS.close();
	PNORMS.close();
	ENERGY.close();
	ENERGYDIFF.close();
	XMAKEMOL.close();
} // end function load()
