#ifndef HARMONIC_OSCILLATOR_3D_H
#define HARMONIC_OSCILLATOR_3D_H

class Harmonic_Oscillator_3d
{
    double **A;
    double **P;

    int n_steps;

    double rho_min;
    double rho_max;
    double h;

    double omega_r;

    double rho;
    double d;
    double e;

    bool interaction;

public:
    Harmonic_Oscillator_3d(int n_step, double rho_max, bool interaction = false,
                           double omega_r = 0);
    void solve (bool eigenvectors = false);
};

#endif // HARMONIC_OSCILLATOR_3D_H
