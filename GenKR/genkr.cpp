#include "hamiltonian.hpp"


//-------- Apply Generalized GenKR boundary conditions ---------
void hamiltonian::deform_box(double dt)
{
    lam += dlam * dt;
    double ipart;
    double tlam = lam(0, 0);
    // For this coordinate system, the pbcs are simple.
    if (fabs(round(lam(0, 0))) + fabs(round(lam(1, 0))) > 0.5)
        reset_count++;
    lam(0, 0) = lam(0, 0) - round(lam(0, 0));
    lam(1, 0) = lam(1, 0) - round(lam(1, 0));

    nPer += lam(0, 0) - tlam;
    alpha_t = omega * lam;
    expAt = alpha_t.exp_diag();
    expmAt = (alpha_t * -1).exp_diag();

    L = S * expAt * Vinv;
    Linv = V * expmAt * Sinv;
}

hamiltonian::hamiltonian(const input &I)
{
    // Simulation box size
    dim = I.dim;
    Ndim_x = I.Ndim_x;
    Ndim_y = I.Ndim_y;
    Ndim_z = I.Ndim_z;
    np = I.Ndim_x * I.Ndim_y * I.Ndim_z;
    a = I.a;
    volume = a * a * a * np;
    dcut = I.dcut;

    // L-J parameters
    eps = I.eps;
    sig = I.sig;
    p1 = 12.;
    p2 = 6.;
	D = -4*eps*(pow(sig/dcut,p1) - pow(sig/dcut,p2));
	C = 4.*eps*(p1*pow(sig/dcut,p1) - p2*pow(sig/dcut,p2))/dcut;


    nChainLength = 0;
    nChainFreq = np; // This must be at least chainLength, =np means a single chain.

    // Track minimum particle distance
    mindist = 2 * a * max(max(Ndim_x, Ndim_y), Ndim_z);
    A = I.A;
    Matrix diag1(3, 3);
    diag1.zeros();
    diag1(0, 0) = Ndim_x * a;
    diag1(1, 1) = Ndim_y * a;
    diag1(2, 2) = Ndim_z * a;
    S = I.S;
    S = S * diag1;
    diag1(0, 0) = 1. / (Ndim_x * a);
    diag1(1, 1) = 1. / (Ndim_y * a);
    diag1(2, 2) = 1. / (Ndim_z * a);
    Sinv = I.Sinv;
    Sinv = Sinv * diag1;

    // Generalized KR boundary conditions
    double tmp1, tmp2, tmp3;
    V.set_size(3, 3);
    tmp1 = 0.73697622909957805;
    tmp2 = 0.59100904850610347;
    tmp3 = 0.32798527760568191;
    V(0, 0) = tmp2;
    V(0, 1) = tmp3;
    V(0, 2) = -tmp1;
    V(1, 0) = -tmp1;
    V(1, 1) = tmp2;
    V(1, 2) = -tmp3;
    V(2, 0) = tmp3;
    V(2, 1) = tmp1;
    V(2, 2) = tmp2;

    Vinv = V.transpose();

    omega.set_size(3, 2);
    omega(0, 0) = 1.619173832;
    omega(0, 1) = -1.177725212;
    omega(1, 0) = -0.441448621;
    omega(1, 1) = 1.619173832;
    omega(2, 0) = -1.177725212;
    omega(2, 1) = -0.441448621;

    Matrix omega_ul = omega.subslice(0, 1, 0, 1);

    lam.set_size(2, 1);
    lam.zeros();
    reset_count = 0;

    // This is our definition, but it is trivial in the constructor.
    alpha_t = omega * lam;

    Matrix tmp(2, 1);
    tmp(0, 0) = A(0, 0);
    tmp(1, 0) = A(1, 1);
    dlam = omega_ul.inverse() * tmp; // FIXME: I've used inverse !

    expAt = alpha_t.exp_diag(); // Should just return identity
    expmAt = (alpha_t * (-1.)).exp_diag();
    L = S * expAt * Vinv;
    Linv = V * expmAt * Sinv;
    initialize_clist(); // note: added to setup clist
}
