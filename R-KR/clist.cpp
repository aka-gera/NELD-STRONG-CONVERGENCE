#include "hamiltonian.hpp"

void hamiltonian::initialize_clist()
{
	region[0] = sqrt(L(0, 0) * L(0, 0) + L(0, 1) * L(0, 1) + L(0, 2) * L(0, 2));
	region[1] = sqrt(L(1, 0) * L(1, 0) + L(1, 1) * L(1, 1) + L(1, 2) * L(1, 2));
	region[2] = sqrt(L(2, 0) * L(2, 0) + L(2, 1) * L(2, 1) + L(2, 2) * L(2, 2)); 
	
	lcmax = static_cast<int>(std::round(a / dcut * Ndim_x));

	cout << "lcmax " << lcmax << " a " << a << endl;
	lc[0] = lcmax;
	lc[1] = lcmax;
	lc[2] = lcmax;
	head = new int[lc[0] * lc[1] * lc[2]];
	lscl = new int[np];
	clear_clist();
}

void hamiltonian::deconstruct_clist()
{
	delete[] head;
	delete[] lscl;
}

void hamiltonian::clear_clist()
{
	for (int i = 0; i < np; i++)
	{
		lscl[i] = -1;
	}

	for (int i = 0; i < lcmax * lcmax * lcmax; i++)
	{
		head[i] = -1;
	}
}

void hamiltonian::add_particles(particle &X)
{
	int c;
	int mc[3]; // holds x,y,z index of a cell
	double da[3];
	Matrix cross = L.crosses();
	min_box_height = min(min(volume / cross.column_norm(0), volume / cross.column_norm(1)), volume / cross.column_norm(2));

	for (int i = 0; i < 3; i++)
	{
		lc[i] = (int)floor(max(min(((volume / cross.column_norm(i)) / dcut), (double)lcmax), 1.));
		// FIXME: add special handling for planar case here.
		if (!(lc[i] >= 3 && lc[i] <= lcmax))
		{
			cout << "Invalid number of boxes " << lc[i] << " of maximum " << lcmax << " in direction " << i << ".  Cell list invalid." << endl;
			cout << "dcut: " << dcut << " box height " << volume / cross.column_norm(i) << endl;
			exit(-1);
		}
	}
	nCells = lc[0] * lc[1] * lc[2];
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
			mc[a] = (int)floor(da[a] * lc[a]);
			mc[a] = (mc[a] + lc[a]) % lc[a];

			if (!(mc[a] >= 0 && mc[a] < lc[a]))
			{
				cout << endl;
				cout << "Bad index in add_particles for particle " << i
					 << " position values:"
					 << ": " << X.q(0, i) << ", " << X.q(1, i) << ", " << X.q(2, i)
					 << "   Particle da vals " << i << ": " << setw(12) << da[0] << ", " << setw(12) << da[1] << ", " << setw(12) << da[2] << endl
					 << "index value for " << a << " :" << mc[a] << ", lc " << a << " value: " << lc[a] << endl;
				cout << endl;
			}
		}
		// translate to cell index, scheme obtained from LL cell MD paper
		c = mc[0] * lc[1] * lc[2] + mc[1] * lc[2] + mc[2];
		if (c > lcmax * lcmax * lcmax)
		{
			cout << endl;
			cout << "trying to access array out of bounds here, c=" << c << ", lcmax^3=" << lcmax * lcmax * lcmax << endl;
			exit(-1);
		}
		// particle linked to previous occupant of head, then head set to particle
		lscl[i] = head[c];
		head[c] = i;
	}
}
