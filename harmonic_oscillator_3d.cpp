#include <cmath>
#include <cstring>
#include "jacobi_algorithm.h"
#include "harmonic_oscillator_3d.h"
#include "lib.h"

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
    rho = rho_min;

    omega_r2 = omega_r*omega_r;

    h = ( rho_max - rho_min )/(n_steps+1);

    e = -1/(h*h);

    if( interact ){
        for( i = 0; i < (n_steps-1); i++ ){
            rho += h;

            A[i][i] = -2*e + omega_r2*rho*rho + 1/rho;
            A[i+1][i] = e;
            A[i][i+1] = e;

            for( j = 0; j < n; j++ ){
                if ( fabs(i-j) > 1 ){
                    A[i][j] = 0;
                }
            }
        }
        A[n_steps-1][n_steps-1] = -2*e + omega_r2*rho*rho + 1/rho;
    }

    else{
        for( i = 0; i < (n_steps-1); i++ ){
            rho += h;

            A[i][i] = -2*e + rho*rho;
            A[i+1][i] = e;
            A[i][i+1] = e;

            for( j = 0; j < n; j++ ){
                if ( fabs(i-j) > 1 ){
                    A[i][j] = 0;
                }
            }
        }
        A[n_steps-1][n_steps-1] = -2*e + rho*rho;
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


void Harmonic_Oscillator_3d::solve(double *eigenvalues, int *counter, const char *method, bool normalize)
{
    int i;
    int j;

    double e_value;
    double *e_list = new double[n+1];

    if( strcmp(method, "jacobi") == 0 ){
        if( eigenvectors ){
            jacobi_algorithm(n, A, P, counter);
        }
        else{
            jacobi_algorithm(n, A, 0, counter);
        }

        for( i = 0; i < n; i++ ){
            eigenvalues[i] = A[i][i];
        }
    }
    else if( strcmp(method, "householder") == 0 ){
        e_value = A[0][1];
        e_list[0] = 0;
        for( i = 0; i < n; i++ ){
            e_list[i+1] = e_value;
            eigenvalues[i] = A[i][i];
        }

        for( i = 0; i < n; i++ ){
            for( j = 0; j < n; j++ ){
                if( i == j ){
                    P[i][j] = 1;
                }
                else{
                    P[i][j] = 0;
                }
            }
        }

        tqli(eigenvalues,e_list,n,P);
    }

    if( normalize ){
        for( i = 0; i < n; i++){
            for( j = 0; j < n; j++){
                P[i][j] = -P[i][j]/sqrt(h);
            }
        }
    }

    delete[] e_list;
}
