#include "harmonic_oscillator_3d.h"
#include "jacobi_algorithm.h"
#include <cmath>

#include <iostream>

using namespace std;

Harmonic_Oscillator_3d::Harmonic_Oscillator_3d(double **A_matrix, double **P_matrix,
                                               bool interaction)
{
    A = A_matrix;
    if( P_matrix != 0 ){
        P = P_matrix;
        eigenvectors = true;
    }
    else{
        eigenvectors = false;
    }
    interact = interaction;
}

void Harmonic_Oscillator_3d::initialize(int n_steps, double rho_max,
                                        double omega_r, double **P_matrix)
{
    n = n_steps;

    if( P_matrix != 0){
        P = P_matrix;
        eigenvectors = true;
    }

    if( omega_r != 0){
        interact = true;
    }

    int i;
    int j;

    double omega_r2;

    rho_min = 0;

    omega_r2 = omega_r*omega_r;

    h = ( rho_max - rho_min )/n_steps;

    e = -1/(h*h);

    if( interact ){
        for( i = 0; i < n_steps; i++ ){
            rho = rho_min + i*h;
            if( rho > 0 ){
                d = -2*e + omega_r2*rho*rho + 1/rho;
            }
            else{
                d = -2*e;
            }

            for( j = 0; j < n_steps; j++){
                if( i == j ){
                        A[i][j] = d;
                }
                else if( fabs( i - j ) == 1){
                    A[i][j] = e;
                }
                else if( fabs( i - j ) > 1){
                    A[i][j] = 0;
                }
            }
        }
    }

    else{
        for( i = 0; i < n_steps; i++ ){
            rho = rho_min + i*h;
            d = -2*e + rho*rho;

            for( j = 0; j < n_steps; j++){
                if( i == j ){
                        A[i][j] = d;
                }
                else if( fabs( i - j ) == 1){
                    A[i][j] = e;
                }
                else if( fabs( i - j ) > 1){
                    A[i][j] = 0;
                }
            }
        }
    }

    if( eigenvectors ){
        for( i = 0; i < n_steps; i++ ){
            for( j = 0; j < n_steps; j++ ){
                if( i == j ){
                    P[i][j] = 1;
                }
                else{
                    P[i][j] = 0;
                }
            }
        }
    }
}


void Harmonic_Oscillator_3d::solve(double *eigenvalues)
{
    int i;

    if( eigenvectors ){
        jacobi_algorithm(n, A, P);
    }
    else{
        jacobi_algorithm(n, A);
    }

    for( i = 0; i < n; i++ ){
        eigenvalues[i] = A[i][i];
    }
}
