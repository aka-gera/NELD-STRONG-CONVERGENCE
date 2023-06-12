#include "algorithm.hpp"
#include "hamiltonian.hpp"
#include <time.h>

int main(void)
{
	//---------- reading data from "input_file" --------------
	time_t start = time(NULL);
	input I;
	I.load();

	//----------- creating structures ------------------------
	algorithm *A;
	A = new algorithm(I);
	A->H = new hamiltonian(I);
	A->Initialize();

	//----------- effective computation ------------------------
	A->run();
	time_t end_t = time(NULL);
	A->H->deconstruct_clist();
	cout << "Total wall time spent running: " << end_t - start << endl;
	delete A->H;
	delete A;

	return EXIT_SUCCESS;
}
