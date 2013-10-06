#ifndef HARMONIC_OSCILLATOR_3D_H
#define HARMONIC_OSCILLATOR_3D_H


class Harmonic_Oscillator_3d
{
    int n;

    double **A;
    double **P;

    double rho_min;
    double h;


    double rho;
    double d;
    double e;

    bool interact;
    bool eigenvectors;

public:
    Harmonic_Oscillator_3d(double **A_matrix, double **P_matrix = 0,
                           bool interaction = false);
    void initialize(int n_steps, double rho_max, double omega_r = 0,
                    double **P_matrix = 0);
    void solve (double *eigenvalues, int *counter = 0, const char *method = "jacobi",
                bool normalize = false);
};

#endif // HARMONIC_OSCILLATOR_3D_H
