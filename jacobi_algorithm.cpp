#include <cmath>
#include "jacobi_algorithm.h"

using namespace std;

void jacobi_algorithm (int n, double **A, double **P, double epsilon)
{
    int k;
    int l;

    double s;
    double c;

    double max_a_ij;

    max_a_ij = max_nondiag(n, A, &k, &l);

    if( P == 0 ){
        while( max_a_ij > epsilon ){
            similarity_transformation(n, A, &k, &l, &c, &s);
            max_a_ij = max_nondiag(n, A, &k, &l);
        }
    }
    else{
        while( max_a_ij > epsilon ){
            similarity_transformation(n, A, &k, &l, &c, &s);
            eigenvectors(n, P, &k, &l, &c, &s);
            max_a_ij = max_nondiag(n, A, &k, &l);
        }
    }
}


void similarity_transformation (int n, double **A, int *k, int *l, double *c, double *s)
{
    int i;

    static double tau;
    static double t;
    static double c2;
    static double s2;
    static double cs;

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

    *c = 1 / sqrt(1 + t*t);
    *s = t*(*c);

    c2 = (*c)*(*c);
    s2 = (*s)*(*s);
    cs = (*c)*(*s);

    a_kk = A[*k][*k];
    a_ll = A[*l][*l];

    A[*k][*k] = a_kk*c2 - 2*A[*k][*l]*cs + a_ll*s2;
    A[*l][*l] = a_ll*c2 + 2*A[*k][*l]*cs + a_kk*s2;


    A[*k][*l] = 0;
    A[*l][*k] = 0;

    for( i = 0; i < n; i ++){
        if( i != *k && i != *l ){
            a_ik = A[i][*k];
            a_il = A[i][*l];

            A[i][*k] = a_ik*(*c) - a_il*(*s);
            A[*k][i] = A[i][*k];
            A[i][*l] = a_il*(*c) + a_ik*(*s);
            A[*l][i] = A[i][*l];
        }
    }
}

void eigenvectors (int n, double **P, int *k, int *l, double *c, double *s)
{
    int i;

    static double *P_k = new double[n];
    static double *P_l = new double[n];

    for( i = 0; i < n; i++ ){
        P_k[i] = P[i][*k];
        P_l[i] = P[i][*l];
    }

    for( i = 0; i < n; i++ ){
        P[i][*k] = *c * P_k[i] - *s * P_l[i];
        P[i][*l] = *s * P_k[i] + *c * P_l[i];
    }
}


double max_nondiag (int n, double **A, int *k, int *l)
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
