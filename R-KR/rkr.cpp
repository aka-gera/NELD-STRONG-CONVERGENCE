#include "hamiltonian.hpp"

//-------- Apply Generalized R-KR boundary conditions ---------
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

    offdiag(0, 0) = cos(rkr_beta * nPer);
    offdiag(0, 1) = -sin(rkr_beta * nPer);
    offdiag(1, 0) = sin(rkr_beta * nPer);
    offdiag(1, 1) = cos(rkr_beta * nPer);

    offdiaginv(0, 0) = cos(rkr_beta * nPer);
    offdiaginv(0, 1) = sin(rkr_beta * nPer);
    offdiaginv(1, 0) = -sin(rkr_beta * nPer);
    offdiaginv(1, 1) = cos(-rkr_beta * nPer);

    L = S * expAt * offdiag * Vinv;
    Linv = V * offdiaginv * expmAt * Sinv;
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

    offdiag.set_size(3, 3);
    offdiag.identity();

    offdiaginv.set_size(3, 3);
    offdiaginv.identity();

    // Rotating KR boundary conditions
    double tmp1, tmp2, tmp3;
    V.set_size(3, 3);
    V(0, 0) = 1.156496717060367;
    V(0, 1) = 0;
    V(0, 2) = -0.311800831557343;
    V(1, 0) = -0.390480732768135;
    V(1, 1) = 0.650274409770356;
    V(1, 2) = 0.724848992205197;
    V(2, 0) = -0.532222319900072;
    V(2, 1) = -0.211155777989163;
    V(2, 2) = 1.272021307822323;

    Vinv.set_size(3, 3);
    Vinv(0, 0) = 0.980218958033226;
    Vinv(0, 1) = 0.065838547165159;
    Vinv(0, 2) = 0.202756101706857;
    Vinv(1, 0) = 0.110919000166463;
    Vinv(1, 1) = 1.305141104609129;
    Vinv(1, 2) = -0.716533262665601;
    Vinv(2, 0) = 0.428542817857063;
    Vinv(2, 1) = 0.244200964032794;
    Vinv(2, 2) = 0.752040220087784;

    omega.set_size(3, 2);
    omega(0, 0) = rkr_mu;
    omega(0, 1) = 0.;
    omega(1, 0) = rkr_mu;
    omega(1, 1) = 0.;
    omega(2, 0) = -2 * rkr_mu;
    omega(2, 1) = 0.;

    lam.set_size(2, 1);
    lam.zeros();
    reset_count = 0;

    alpha_t = omega * lam;
    dlam.set_size(2, 1);
    dlam.zeros();
    dlam(0, 0) = A(0, 0) / omega(0, 0);

    expAt = alpha_t.exp_diag(); // Should just return identity
    expmAt = (alpha_t * (-1.)).exp_diag();
    L = S * expAt * Vinv;
    Linv = V * expmAt * Sinv;
    initialize_clist(); // note: added to setup clist
}
