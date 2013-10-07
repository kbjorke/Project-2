#include <iostream>
#include <iomanip>
#include <fstream>
#include "time.h"
#include "output_functions.h"
#include "harmonic_oscillator_3d.h"

using namespace std;

double getUnixTime();

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

void output_eigenvalues(int n, double omega_r, double rho_max,
                        int m, const char *method)
{
    int j;
    int k;

    int counter;

    double *index_list = new double[m];

    double *eigenvalues = new double[n];

    double **A = new double*[n];
    for( j = 0; j < n; j++ ){
         A[j] = new double[n];
    }

    Harmonic_Oscillator_3d harmonic_oscillator (A);

    harmonic_oscillator.initialize(n, rho_max, omega_r);

    harmonic_oscillator.solve(eigenvalues,&counter, method);


    eigenstates_index(m, index_list, n, eigenvalues);


    cout << "n : " << n << "  omega_r : " << omega_r << "  rho_max : " <<
            rho_max << endl;

    for( j = 0; j < m; j++ ){
        k = index_list[j];
        cout << "lambda_" << j << ": ";
        cout << setprecision(14) << scientific << eigenvalues[k] << endl;
    }

    for( j = 0; j < n; j++ ){
        delete[] A[j];
    }
    delete[] A;

    delete[] eigenvalues;

    delete[] index_list;
}

void output_stability(int n_start, int n_stop, double resolution, double rho_max,
                      int m, const char *method)
{
    int n;
    int i;
    int j;
    int k;


    double start,end,time;

    double h;

    int counter;

    h = ( n_stop - n_start ) / (resolution-1);

    double *index_list = new double[m];

    // Write to file:
    fstream myfile;
    myfile.open("output_stability.txt", ios::out);

    myfile << resolution << '\t' << rho_max << '\t' << m << endl;

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

        start = getUnixTime();
        harmonic_oscillator.solve(eigenvalues,&counter, method);
        end = getUnixTime();

        time = end-start;
        cout << n << ":  " << time << " s" << endl;

        eigenstates_index(m, index_list, n, eigenvalues);


        for( j = 0; j < m; j++ ){
            k = index_list[j];
            myfile << setprecision(14) << scientific << '\t' << eigenvalues[k];
        }
        myfile << '\t' << counter << endl;

        for( j = 0; j < n; j++ ){
            delete[] A[j];
        }
        delete[] A;

        delete[] eigenvalues;
    }

    myfile.close();

    delete[] index_list;
}

void output_eigenstates(int m, int n_steps, double rho_max,
                        double omega_r, const char *method)
{
    int i;
    int j;
    int k;

    int counter;

    double *index_list = new double[m];

    // Write to file:
    fstream myfile;
    myfile.open("output_eigenstates.txt", ios::out);

    myfile << m << '\t' << n_steps << '\t' << rho_max << '\t' << omega_r
           << '\t' << endl;

    double *eigenvalues = new double[n_steps];

    double **A = new double*[n_steps];
    for( j = 0; j < n_steps; j++ ){
         A[j] = new double[n_steps];
    }

    double **P = new double*[n_steps];
    for( j = 0; j < n_steps; j++ ){
         P[j] = new double[n_steps];
    }

    Harmonic_Oscillator_3d harmonic_oscillator (A, P);

    harmonic_oscillator.initialize(n_steps, rho_max, omega_r);

    harmonic_oscillator.solve(eigenvalues,&counter, method,true);

    eigenstates_index(m, index_list, n_steps, eigenvalues);

    for( i = 0; i < n_steps; i++ ){
        for( j = 0; j < m; j++ ){
                k = index_list[j];
                myfile << setprecision(14) << scientific << P[i][k] << '\t';
        }
            myfile << endl;
    }

    for( j = 0; j < n_steps; j++ ){
        delete[] A[j];
    }
    delete[] A;
    for( j = 0; j < n_steps; j++ ){
        delete[] P[j];
    }
    delete[] P;

    delete[] eigenvalues;
    delete[] index_list;

    myfile.close();
}


double getUnixTime()
{
    /* Function to get accurate time, prescission to about some microseconds.
     *
     * Important: to get the function clock_gettime() to work link to library
     * -lrt when linking.
     * */
    struct timespec tv;

    // Get time:
    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    // Returns time in seconds:
    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}

