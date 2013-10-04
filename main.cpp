#include <unittest++/UnitTest++.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "jacobi_algorithm.h"
#include "harmonic_oscillator_3d.h"

#include <algorithm>

using namespace std;

void eigenstates_index(int m, double *index_list, int n, double *eigenvalues);
void output_stability(int n_start, int n_stop, int resolution, double rho_max,
                      int m = 3);

int main()
{

    output_stability(140, 141, 2, 4.5);


    /*
    int n;
    int i;
    int j;

    double rho_max;

    double *index_list = new double[3];

    rho_max = 4.5;
    n = 140;

    double *eigenvalues = new double[n];
    int counter;

    double **A = new double*[n];
    for( i = 0; i < n; i++ ){
        A[i] = new double[n];
    }

    double **P = new double*[n];
    for( i = 0; i < n; i++ ){
        P[i] = new double[n];
    }


    Harmonic_Oscillator_3d harmonic_oscillator (A);

    harmonic_oscillator.initialize(n, rho_max);

    harmonic_oscillator.solve(eigenvalues,&counter);

    eigenstates_index(3, index_list, n, eigenvalues);

    for( i = 0; i < 3; i++ ){
        j = index_list[i];
        cout << eigenvalues[j] << " " << index_list[i] << endl;
    }

    sort(eigenvalues, eigenvalues+n);

    for( i = 0; i < n; i++ ){
        cout << eigenvalues[n-1-i] << endl;
    }

    cout << endl << counter << endl;




    for( i = 0; i < n; i++ ){
        delete[] A[i];
    }
    delete[] A;



    for( i = 0; i < n; i++ ){
        delete[] P[i];
    }
    delete[] P;
    */


    //return UnitTest::RunAllTests();
}

void eigenstates_index(int m, double *index_list, int n, double *eigenvalues)
{
    int i;
    int j;
    int index;

    double min_value;
    double last_value;

    index = 0;
    last_value = 0;

    for( i = 0; i < m; i++ ){
        min_value = 1e10;
        for( j = 0; j < n; j++ ){
            if( ( eigenvalues[j] < min_value ) && ( eigenvalues[j] > last_value ) ){
                min_value = eigenvalues[j];
                index = j;
            }
        }
        index_list[i] = index;
        last_value = eigenvalues[index];
    }
}


void output_stability(int n_start, int n_stop, int resolution, double rho_max, int m)
{
    int n;
    int i;
    int j;
    int k;

    int h;

    int counter;

    h = ( n_stop - n_start ) / (resolution-1);

    double *index_list = new double[m];

    // Write to file:
    fstream myfile;
    myfile.open("output_stability.txt", ios::out);

    for(i=0; i < resolution; i++){
        n = n_start + i*h;
        cout << n << endl;

        myfile << n;

        double *eigenvalues = new double[n];

        double **A = new double*[n];
        for( j = 0; j < n; j++ ){
             A[j] = new double[n];
        }

        Harmonic_Oscillator_3d harmonic_oscillator (A);

        harmonic_oscillator.initialize(n, rho_max);

        harmonic_oscillator.solve(eigenvalues,&counter);

        eigenstates_index(m, index_list, n, eigenvalues);

        for( j = 0; j < m; j++ ){
            k = index_list[j];
            cout << k << "  " << eigenvalues[k] << endl;
            myfile << '\t' << eigenvalues[k];
        }
        myfile << '\t' << counter << endl;

        for( j = 0; j < n; j++ ){
            delete[] A[j];
        }
        delete[] A;

        delete[] eigenvalues;
    }
    myfile.close();
}

void output_eigenstates(){

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
