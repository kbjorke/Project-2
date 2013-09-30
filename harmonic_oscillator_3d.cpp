#include "harmonic_oscillator_3d.h"

Harmonic_Oscillator_3d::Harmonic_Oscillator_3d(int n_step, double rho_max,
                                               bool interaction, double omega_r)
{
    rho_min = 0;
    h = ( rho_max - rho_min )/ n_steps;
}

void Harmonic_Oscillator_3d::solve(bool eigenvectors)
{
    int i;
    int j;

    double **A = new double*[n_steps];
    for( i = 0; i < n_steps; i++ ){
        A[i] = new double[n_steps];
    }

    e = -1/(h*h);

    A[0][0] = -2*e;

    if( interaction ){
        for( i = 1; i < n_steps; i++ ){
            rho = rho_min + i*h;
            d = 2/(h*h) + omega_r*omega_r*rho*rho + 1/rho;

            A[i][i] = d;
            A[i-1][i] = e;
            A[i][i-1] = e;
        }
    }
    else{
        for( i = 1; i < n_steps; i++ ){
            rho = rho_min + i*h;
            d = 2/(h*h) + rho*rho;

            A[i][i] = d;
            A[i-1][i] = e;
            A[i][i-1] = e;
        }
    }
}
