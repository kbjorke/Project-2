#include <cmath>
#include "jacobi_algorithm.h"

#include <iostream>

using namespace std;

void jacobi_algorithm (int n, double **A, double epsilon)
{
    int k;
    int l;

    double max_a_ij;

    max_a_ij = max_nondiag(n, A, &k, &l);

    while( max_a_ij > epsilon ){
        similarity_transformation(n, A, &k, &l);
        max_a_ij = max_nondiag(n, A, &k, &l);
    }
}


void similarity_transformation(int n, double **A, int *k, int *l)
{
    int i;

    static double tau;
    static double t;
    static double s;
    static double c;

    static double a_kk;
    static double a_ll;
    static double a_ik;
    static double a_il;

    tau = ( A[*l][*l] - A[*k][*k] ) / ( 2*A[*k][*l] );

    if( tau > 0 ){
        t = 1 / (tau + sqrt(1 + tau*tau));
    }
    else{
        t = -1 / (-tau + sqrt(1 + tau*tau));
    }

    c = 1 / sqrt(1 + t*t);
    s = t*c;

    a_kk = A[*k][*k];
    a_ll = A[*l][*l];

    A[*k][*k] = a_kk*c*c - 2*A[*k][*l]*c*s + a_ll*s*s;
    A[*l][*l] = a_ll*c*c + 2*A[*k][*l]*c*s + a_kk*s*s;


    A[*k][*l] = 0;
    A[*l][*k] = 0;

    for( i = 0; i < n; i ++){
        if( i != *k && i != *l ){
            a_ik = A[i][*k];
            a_il = A[i][*l];

            A[i][*k] = a_ik*c - a_il*s;
            A[*k][i] = A[i][*k];
            A[i][*l] = a_il*c + a_ik*s;
            A[*l][i] = A[i][*l];
        }
    }
}


double max_nondiag(int n, double **A, int *k, int *l)
{
    int i;
    int j;

    double max_element = 0;

    for( i = 0; i < n; i++ ){
        for( j = 0; j < n; j ++ ){
            if( i != j ){
                if( fabs(A[i][j]) > max_element ){
                    *k = i;
                    *l = j;
                    max_element = fabs(A[i][j]);
                }
            }
        }
    }
    return max_element;
}
