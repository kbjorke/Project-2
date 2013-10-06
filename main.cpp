#include <unittest++/UnitTest++.h>
#include <iostream>
#include <cmath>
#include "jacobi_algorithm.h"
#include "harmonic_oscillator_3d.h"
#include "output_functions.h"

#include <algorithm>

using namespace std;

int main()
{

    //output_stability(50, 100, 10, 6, 3, "householder");
    output_eigenstates(3, 300, 6, 0, "householder");


    //return UnitTest::RunAllTests();
}




TEST(max_nondiag_func)
{
    int i;
    int size = 2;
    double max_element;

    double **matrix = new double*[size];
    for( i = 0; i < size; i++ ){
        matrix[i] = new double[size];
    }

    matrix[0][0] = 1;
    matrix[0][1] = 2;
    matrix[1][0] = 4;
    matrix[1][1] = 3;

    int k;
    int l;

    max_element = max_nondiag(size, matrix, &k, &l);

    CHECK( ( max_element == 4 ) && k == 1 && l == 0);

    for( i = 0; i < size; i++ ){
         delete[] matrix[i];
    }
    delete[] matrix;
}

TEST(jacobi_algorithm)
{
    int i;
    int n = 3;
    int counter;

    double epsilon = 1e-8;

    double **A = new double*[n];
    for( i = 0; i < n; i++ ){
        A[i] = new double[n];
    }

    A[0][0] = 6;
    A[0][1] = -2;
    A[0][2] = -1;
    A[1][0] = -2;
    A[1][1] = 6;
    A[1][2] = -1;
    A[2][0] = -1;
    A[2][1] = -1;
    A[2][2] = 5;

    jacobi_algorithm(n, A, 0, &counter, epsilon);

    /*
    for( i = 0; i < n; i++ ){
        for( j = 0; j < n; j++ ){
            cout << A[i][j] << '\t';
        }
        cout << endl;
    }
    */

    CHECK( ( A[0][0] - 3 < epsilon ) && ( A[0][1] < epsilon ) && ( A[0][2] < epsilon ) &&
            ( A[1][0] < epsilon ) && ( A[1][1] - 8 < epsilon ) && ( A[1][2] < epsilon ) &&
            ( A[2][0] < epsilon ) && ( A[2][1] < epsilon ) && ( A[2][2] - 6 < epsilon ) );


    for( i = 0; i < n; i++ ){
        delete[] A[i];
    }
    delete[] A;
}

TEST(jacobi_eigenvectors)
{
    int i;
    int j;
    int n = 3;

    int counter;

    double epsilon = 1e-8;

    double **A = new double*[n];
    for( i = 0; i < n; i++ ){
        A[i] = new double[n];
    }
    double **P = new double*[n];
    for( i = 0; i < n; i++ ){
        P[i] = new double[n];
    }


    A[0][0] = 6;
    A[0][1] = -2;
    A[0][2] = -1;
    A[1][0] = -2;
    A[1][1] = 6;
    A[1][2] = -1;
    A[2][0] = -1;
    A[2][1] = -1;
    A[2][2] = 5;

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

    jacobi_algorithm(n, A, P, &counter, 1e-8);

    CHECK( ( fabs(P[0][0]) - 1/sqrt(3) < epsilon ) &&
           ( fabs(P[1][0]) - 1/sqrt(3) < epsilon ) &&
           ( fabs(P[2][0]) - 1/sqrt(3) < epsilon ) &&

           ( fabs(P[0][1]) - 1/sqrt(2) < epsilon ) &&
           ( fabs(P[1][1]) - 1/sqrt(2) < epsilon ) &&
           ( fabs(P[2][1]) - 0 < epsilon ) &&

           ( fabs(P[0][2]) - 1/sqrt(6) < epsilon ) &&
           ( fabs(P[1][2]) - 1/sqrt(6) < epsilon ) &&
           ( fabs(P[2][2]) - 2/sqrt(6) < epsilon ) );

    for( i = 0; i < n; i++ ){
        delete[] A[i];
    }
    delete[] A;
    for( i = 0; i < n; i++ ){
        delete[] P[i];
    }
    delete[] P;
}

TEST(Harmonic_ocillator_3d_initialize)
{
    int n;
    int i;
    int j;

    double rho_max;

    rho_max = 4.5;
    n = 8;

    double **A = new double*[n];
    for( i = 0; i < n; i++ ){
        A[i] = new double[n];
    }
    double **P = new double*[n];
    for( i = 0; i < n; i++ ){
        P[i] = new double[n];
    }

    Harmonic_Oscillator_3d harmonic_oscillator (A, P);

    harmonic_oscillator.initialize(n, rho_max);


    for( i = 0; i < n; i++ ){
        for( j = 0; j < n; j++ ){
            cout << A[i][j] << '\t';
        }
        cout << endl;
    }
    cout << endl;
    for( i = 0; i < n; i++ ){
        for( j = 0; j < n; j++ ){
            cout << P[i][j] << '\t';
        }
        cout << endl;
    }


    for( i = 0; i < n; i++ ){
        delete[] A[i];
    }
    delete[] A;
    for( i = 0; i < n; i++ ){
        delete[] P[i];
    }
    delete[] P;
}


TEST(Harmonic_ocillator_3d_solve)
{
}
